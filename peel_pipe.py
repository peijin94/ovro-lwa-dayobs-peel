#!/usr/bin/env python3
"""
Minimal pipeline script for LWA data processing
Complete workflow: raw MS -> CASA applycal -> DP3 flag/avg -> wsclean -> gaincal -> applycal
"""
import subprocess, sys, os
import time
from pathlib import Path
import shutil
import argparse
import wsclean_imaging
import h5py
import numpy as np
import config
from source_list import get_time_mjd, get_Sun_RA_DEC, mask_far_Sun_sources, src_coord_lst

PIPELINE_SCRIPT_DIR = Path(__file__).parent
EXECUTABLE_DIR = Path(__file__).parent / "exe"

def run_casa_applycal(input_ms, gaintable):
    """Apply CASA bandpass calibration"""
    print(f"Step : CASA applycal - {input_ms}")
    start_time = time.time()
    try:
        subprocess.run(["python3", str(EXECUTABLE_DIR / "flagant_applybp.py"), str(input_ms), str(gaintable)], check=True)
        elapsed = time.time() - start_time
        print(f"✓ CASA applycal completed ({elapsed:.1f}s)")
    except subprocess.CalledProcessError as e:
        elapsed = time.time() - start_time
        print(f"✗ CASA applycal failed after {elapsed:.1f}s: {e.stdout}")
        sys.exit(1)

def run_combine_ms(input_ms_list, output_ms, datacolumn="DATA"):
    """DP3 combine multiple MS files"""
    print(f"Step : DP3 combine MS - {len(input_ms_list)} files -> {output_ms}")
    start_time = time.time()
    
    output_path = Path(output_ms)
    
    # Create msin list string: [ms1,ms2,ms3,...]
    msin_str = "[" + ",".join([str(Path(ms)) for ms in input_ms_list]) + "]"
    
    parset_content = f"""msin={msin_str}
        msout={str(output_path)}
        msin.datacolumn={datacolumn}
        steps=[]
        """
    
    cmd = ["DP3", *parset_content.split()]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        elapsed = time.time() - start_time
        print(f"✓ DP3 combine completed ({elapsed:.1f}s): {output_ms}")
    except subprocess.CalledProcessError as e:
        elapsed = time.time() - start_time
        print(f"✗ DP3 combine failed after {elapsed:.1f}s: {e.stdout}")
        sys.exit(1)

def run_dp3_flag_avg(input_ms, output_ms, strategy_file=None):
    """DP3 flagging and frequency averaging"""
    print(f"Step : DP3 flag/avg - {input_ms} -> {output_ms}")
    start_time = time.time()
    
    input_path = Path(input_ms)
    output_path = Path(output_ms)
    
    if strategy_file is not None:
        strategy_file = Path(strategy_file)
        strategy_dest = Path("/data") / strategy_file.name
        shutil.copy(strategy_file, strategy_dest)
        strategy_file_path_str = str(strategy_dest)
    else:
        strategy_file_path_str = "/usr/local/share/linc/rfistrategies/lofar-default.lua"

    print(f"Strategy file: {strategy_file_path_str}")

    # Create DP3 parset - use simple filenames since we're in /data
    parset_content = f"""msin={str(input_path)} 
        msout={str(output_path)}
        msin.datacolumn=CORRECTED_DATA
        steps=[flag,avg]
        flag.type=aoflagger
        flag.strategy={strategy_file_path_str}
        avg.type=averager
        avg.freqstep={config.init_avg_n_freq}
        """
#flag.keepstatistics=false

    cmd = ["DP3", *parset_content.split()]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        elapsed = time.time() - start_time
        print(f"✓ DP3 flag/avg completed ({elapsed:.1f}s): {output_ms}")
    except subprocess.CalledProcessError as e:
        elapsed = time.time() - start_time
        print(f"✗ DP3 flag/avg failed after {elapsed:.1f}s: {e.stdout}")
        sys.exit(1)

def run_wsclean_imaging(input_ms, output_prefix="image", auto_pix_fov=True, **kwargs):
    """WSClean imaging"""
    print(f"Step : WSClean imaging - {input_ms}")
    start_time = time.time()
    
    input_path = Path(input_ms)
    
    # Generate WSClean command using utils
    wsclean_cmd_str = wsclean_imaging.make_wsclean_cmd(
        msfile=input_path,  imagename=output_prefix,
        auto_pix_fov=auto_pix_fov, **kwargs)
    
    # Split the command string into a list for subprocess
    wsclean_args = wsclean_cmd_str.split()[1:]  # Remove 'wsclean' from the beginning
    
    cmd = ["wsclean"] + wsclean_args + [str(input_path)]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        elapsed = time.time() - start_time
        print(f"✓ WSClean imaging completed ({elapsed:.1f}s): {output_prefix}*.fits")
    except subprocess.CalledProcessError as e:
        elapsed = time.time() - start_time
        print(f"✗ WSClean imaging failed after {elapsed:.1f}s: {e.stdout, e.stderr}")
        sys.exit(1)

def run_gaincal(input_ms, solution_fname="solution.h5", cal_type="diagonalphase"):
    """DP3 gain calibration"""
    print(f"Step : DP3 gaincal - {input_ms}")
    start_time = time.time()
    
    input_path = Path(input_ms)
    parset_content = f"""msin={str(input_path)}
        showprogress=False
        verbosity="quiet"
        steps=[gaincal]
        msout=.
        gaincal.solint=0
        gaincal.caltype={cal_type}
        gaincal.uvlambdamin=30
        gaincal.maxiter=500
        gaincal.tolerance=1e-5
        gaincal.usemodelcolumn=true
        gaincal.modelcolumn=MODEL_DATA
        gaincal.parmdb={solution_fname}
        """
    
    cmd = ["DP3", *parset_content.split()]
    
    try:
        res = subprocess.run(cmd, check=True, capture_output=True, text=True)
        elapsed = time.time() - start_time
        print(f"✓ DP3 gaincal completed ({elapsed:.1f}s): solution.h5")
    except subprocess.CalledProcessError as e:
        elapsed = time.time() - start_time
        print(f"✗ DP3 gaincal failed after {elapsed:.1f}s: {e, e.stdout, e.stderr}")
        sys.exit(1)
    

def reset_solution_outliers(h5fname, N_sigma=3, reset=True):
    start_time = time.time()
    with h5py.File(h5fname, 'r') as f:
        amp_val = f["sol000"]["amplitude000"]["val"][:]
        weight_val = f["sol000"]["amplitude000"]["weight"][:]
        
    for i in range(amp_val.shape[0]): # time
        for j in range(amp_val.shape[1]): # freq
            for k in range(amp_val.shape[3]): # pol
                outliers = np.where(
                    (amp_val[i,j,:,k] > np.nanmean(amp_val[i,j,:,k]) + N_sigma * np.nanstd(amp_val[i,j,:,k]))
                    | (amp_val[i,j,:,k] < np.nanmean(amp_val[i,j,:,k]) - N_sigma * np.nanstd(amp_val[i,j,:,k]))
                )[0]
                if reset:
                    amp_val[i,j,outliers,k] = np.nan
                    weight_val[i,j,outliers,k] = 0
                else:
                    amp_val[i,j,outliers,k] = 1
                    weight_val[i,j,outliers,k] = 0

    with h5py.File(h5fname ,'a') as f:
        f["sol000"]["amplitude000"]["val"][:] = amp_val
        f["sol000"]["amplitude000"]["weight"][:] = weight_val

    elapsed = time.time() - start_time
    print(f"✓ Reset solution outliers completed ({elapsed:.1f}s): {h5fname}")
    return h5fname

def run_applycal_dp3(input_ms,  output_ms, solution_fname="solution.h5", cal_entry_lst=["phase"]):
    """Apply DP3 calibration solutions"""
    print(f"Step : DP3 applycal - {input_ms} -> {output_ms}")
    start_time = time.time()
    
    input_path = Path(input_ms)
    output_path = Path(output_ms)
    
    parset_content = f"""msin={str(input_path)}
        msout={str(output_path)}
        steps=[applycal]
        showprogress=False
        verbosity="quiet"
        applycal.parmdb={str(solution_fname)}
        applycal.steps=[{','.join(cal_entry_lst)}] \n
        """
    for cal_entry in cal_entry_lst:
        parset_content += f"applycal.{cal_entry}.correction={cal_entry}000 \n"

    cmd = ["DP3", *parset_content.split()]
    try:
        subprocess.run(cmd, check=True)
        elapsed = time.time() - start_time
        print(f"✓ DP3 applycal completed ({elapsed:.1f}s): {output_ms}")
    except subprocess.CalledProcessError as e:
        elapsed = time.time() - start_time
        print(f"✗ DP3 applycal failed after {elapsed:.1f}s: {e.stdout} {e.stderr}")
        sys.exit(1)
    

def run_dp3_gaincal_Ateam_source(input_ms, skymodel_fname, solution_fname="solution.h5", sources="CygA"):
    """DP3 predict"""
    print(f"Step : DP3 predict - {input_ms} with {skymodel_fname} and {sources}")
    start_time = time.time()
    
    input_path = Path(input_ms)
    source_path = Path(skymodel_fname)
    
    parset_content = f"""msin={str(input_path)}
        showprogress=False
        verbosity="quiet"
        msout=.
        steps=[cal0]
        cal0.type=gaincal
        cal0.elementmodel=lwa
        cal0.solint=0
        cal0.caltype=diagonalphase
        cal0.uvlambdamin=30
        cal0.maxiter=500
        cal0.tolerance=1e-5
        cal0.modelcolumn=MODEL_DATA
        cal0.parmdb={str(solution_fname)}
        cal0.sourcedb={str(source_path)}
        cal0.sources={sources}
        """

    cmd = ["DP3", *parset_content.split()]

    try:
        print(f"Running command: {cmd}")
        res = subprocess.run(cmd, check=True)
        elapsed = time.time() - start_time
        print(f"Output: {res.stdout}")
        print(f"✓ DP3 pred-gaincal completed ({elapsed:.1f}s): {input_ms} with {skymodel_fname} and {sources}")
    except subprocess.CalledProcessError as e:
        elapsed = time.time() - start_time
        print(f"✗ DP3 pred-gaincal failed after {elapsed:.1f}s: {e, e.stdout, e.stderr}")
        sys.exit(1)

def run_dp3_subtract(input_ms, output_ms, skymodel_fname, sources=[]):
    """DP3 subtract"""
    print(f"Step : DP3 subtract - {input_ms} -> {output_ms}")
    start_time = time.time()
    
    input_path = Path(input_ms)
    output_path = Path(output_ms)
    source_path = Path(skymodel_fname)
    
    parset_content = f"""msin={str(input_path)}
        showprogress=False
        verbosity="quiet"
        msout={str(output_path)}
        steps=[predict]
        predict.elementmodel=lwa
        predict.type=predict
        predict.sourcedb={str(source_path)}
        predict.sources={','.join(sources)}
        predict.operation=subtract
        """
        
    cmd = ["DP3", *parset_content.split()]
    
    try:
        subprocess.run(cmd, check=True)
        elapsed = time.time() - start_time
        print(f"✓ DP3 subtract completed ({elapsed:.1f}s): {output_ms}")
    except subprocess.CalledProcessError as e:
        elapsed = time.time() - start_time
        print(f"✗ DP3 subtract failed after {elapsed:.1f}s: {e, e.stdout, e.stderr}")
        sys.exit(1)
    

def phaseshift_to_RA_DEC(ms_file, output_ms, ra_deg, dec_deg):
    """Phase shift MS to Sun's coordinates using DP3 PhaseShift step."""
    ms_path = Path(ms_file)
    output_path = Path(output_ms)
    if not ms_path.exists():
        raise FileNotFoundError(f"MS file not found: {ms_path}")
    start_time = time.time()

    # Create parset content - use simple filenames
    parset_content = f"""msin={str(ms_path)}
        msout={str(output_path)}
        showprogress=False
        verbosity="quiet"
        steps=[phaseshift]
        phaseshift.type=phaseshift
        phaseshift.phasecenter=[{ra_deg},{dec_deg}]
        """
    cmd = ["DP3", *parset_content.split()]
    
    try:
        subprocess.run(cmd, check=True)
        elapsed = time.time() - start_time  
        print(f"✓ DP3 phase shift completed ({elapsed:.1f}s): {output_path}")
        return str(output_path)

    except subprocess.CalledProcessError as e:
        print(f"✗ DP3 phase shift failed after {elapsed:.1f}s: {e.stdout}")
        sys.exit(1)
    
def phaseshift_to_sun(ms_file, output_ms):
    """Phase shift MS to Sun's coordinates using DP3 PhaseShift step."""
    time_mjd = get_time_mjd(str(ms_file))
    sun_ra, sun_dec = get_Sun_RA_DEC(time_mjd)
    return phaseshift_to_RA_DEC(ms_file, output_ms, sun_ra, sun_dec)


def run_dp3_avg(input_ms, output_ms, freq_step=4):
    """DP3 frequency averaging"""
    print(f"Step : DP3 frequency averaging - {input_ms} -> {output_ms}")
    start_time = time.time()
    
    input_path = Path(input_ms)
    output_path = Path(output_ms)
    
    parset_content = f"""msin={str(input_path)}
        msout={str(output_path)}
        steps=[avg]        
        showprogress=False
        verbosity="quiet"
        avg.type=averager
        avg.freqstep={freq_step}
        """
        
    cmd = ["DP3", *parset_content.split()]
    
    try:
        subprocess.run(cmd, check=True)
        elapsed = time.time() - start_time
        print(f"✓ DP3 frequency averaging completed ({elapsed:.1f}s): {output_ms}")
    except subprocess.CalledProcessError as e:
        elapsed = time.time() - start_time
        print(f"✗ DP3 frequency averaging failed after {elapsed:.1f}s: {e}")
        sys.exit(1)


def run_reset_zenith(input_ms):
    """Reset zenith with https://wsclean.readthedocs.io/en/latest/chgcentre.html"""
    print(f"Step : Reset zenith - {input_ms}")
    start_time = time.time()
    
    input_path = Path(input_ms)
    cmd = ["chgcentre", "-zenith", str(input_path)]
    try:
        subprocess.run(cmd, check=True)
        elapsed = time.time() - start_time
        print(f"✓ Reset zenith completed ({elapsed:.1f}s): {input_ms}")
    except subprocess.CalledProcessError as e:
        elapsed = time.time() - start_time
        print(f"✗ Reset zenith failed after {elapsed:.1f}s: {e}")
        sys.exit(1)



def run_wsclean_predict(target_ms, image_name):
    """WSClean predict"""
    print(f"Step : WSClean predict - {target_ms} -> {image_name}")
    start_time = time.time()
    cmd = ["wsclean", "-predict", "-name", str(image_name), str(target_ms)]
    try:
        subprocess.run(cmd, check=True)
        elapsed = time.time() - start_time
        print(f"✓ WSClean predict completed ({elapsed:.1f}s): {image_name}")
    except subprocess.CalledProcessError as e:
        elapsed = time.time() - start_time
        print(f"✗ WSClean predict failed after {elapsed:.1f}s: {e}")
        sys.exit(1)


def run_casa_uvsub(target_ms):
    """CASA uvsub"""
    from casatools import uvsub
    uvsub(vis=target_ms, reverse=False)


def run_peel_pipeline(ms_list, bcal_list, output_prefix="peel", output_dir=None):
    """Run peeling pipeline for multiple wideband MS files
    
    Args:
        ms_list: List of input measurement set paths
        bcal_list: List of corresponding bandpass calibration tables
        output_prefix: Prefix for output files
        output_dir: Directory for output files (default: parent of first MS)
    """
    
    pipeline_start = time.time()
    
    # Validate inputs
    if len(ms_list) != len(bcal_list):
        print(f"Error: Number of MS files ({len(ms_list)}) must match number of bandpass tables ({len(bcal_list)})")
        sys.exit(1)
    
    for ms in ms_list:
        if not Path(ms).exists():
            print(f"Error: MS file not found: {ms}")
            sys.exit(1)
    
    for bcal in bcal_list:
        if not Path(bcal).exists():
            print(f"Error: Bandpass table not found: {bcal}")
            sys.exit(1)
    
    # Set output directory
    if output_dir is None:
        output_dir = Path(ms_list[0]).parent
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
    
    # Extract timestamp from first MS filename
    first_ms_name = Path(ms_list[0]).stem
    timestamp = first_ms_name.split('_')[0] + '_' + first_ms_name.split('_')[1]
    
    print("="*60)
    print("LWA Peeling Pipeline - Wideband Processing")
    print("="*60)
    print(f"Input MS files: {len(ms_list)}")
    for i, (ms, bcal) in enumerate(zip(ms_list, bcal_list)):
        print(f"  [{i+1}] {Path(ms).name} -> {Path(bcal).name}")
    print(f"Output directory: {output_dir}")
    print(f"Output prefix: {output_prefix}")
    print("="*60)
    
    # Step 1: Apply bandpass calibration to all MS files
    applied_bp_ms_list = []
    for i, (ms, bcal) in enumerate(zip(ms_list, bcal_list)):
        run_casa_applycal(ms, bcal)
        flag_ms = ms.replace(".ms", "_flag.ms")
        run_dp3_flag_avg(ms, flag_ms, strategy_file=PIPELINE_SCRIPT_DIR / "lua" / "LWA_sun_PZ.lua")
        applied_bp_ms_list.append(str(flag_ms))
    
    # Step 2: Combine MS files
    combined_ms = output_dir / f"{timestamp}_combined.ms"
    print(f"\nCombining {len(applied_bp_ms_list)} MS files...")
    run_combine_ms(applied_bp_ms_list, str(combined_ms))
    
    # Step 3: Initial imaging
    #print(f"\nGenerating initial image...")
    #run_wsclean_imaging(combined_ms, str(output_dir / f"{output_prefix}_initial"), 
    #                   niter=8000, mgain=0.8, horizon_mask=3, auto_pix_fov=False)


    # step 4: phase to coordinate of cyg-A
    src_name = "CygA"
    cyg_a_ms = output_dir / f"{output_prefix}_{src_name}.ms"
    phaseshift_to_RA_DEC(combined_ms, cyg_a_ms, src_coord_lst[src_name][0], src_coord_lst[src_name][1])

    # step 5: frequency averaging [for decorrelation]
    cyg_a_avg_ms = output_dir / f"{output_prefix}_{src_name}_avg.ms"
    run_dp3_avg(cyg_a_ms, cyg_a_avg_ms, freq_step=96)

    # step 6: cal
    skymodel_fname =  PIPELINE_SCRIPT_DIR / "skymodel" / "LOFAR-A-team.skymodel"
    run_dp3_gaincal_Ateam_source(cyg_a_avg_ms, skymodel_fname, solution_fname=f"sol_{src_name}.h5", sources=src_name)

    # step 7: applycal
    applycal_ms = output_dir / f"{output_prefix}_{src_name}_applycal.ms"
    applycal_avg_ms = output_dir / f"{output_prefix}_{src_name}_applycal_avg.ms"
    run_applycal_dp3(cyg_a_ms, applycal_ms, solution_fname=f"sol_{src_name}.h5")
    run_applycal_dp3(cyg_a_avg_ms, applycal_avg_ms, solution_fname=f"sol_{src_name}.h5")

    # step 8: selfcal
    run_wsclean_imaging(applycal_avg_ms, str(output_dir / f"{output_prefix}_{src_name}_avg"), 
            niter=500, mgain=0.9, horizon_mask=3, auto_pix_fov=False, size=128, scale='1arcmin', save_source_list=True, no_negative=True,
            join_channels=True, channels_out=4, fit_spectral_pol=2)
    run_gaincal(applycal_avg_ms, solution_fname=f"sol_self_{src_name}.h5")

    applycal_avg_self_ms = output_dir / f"{output_prefix}_{src_name}_applycal_avg_self.ms"
    applycal_self_ms = output_dir / f"{output_prefix}_{src_name}_applycal_self.ms"
    run_applycal_dp3(applycal_avg_ms, applycal_avg_self_ms, solution_fname=f"sol_self_{src_name}.h5")
    run_applycal_dp3(applycal_ms, applycal_self_ms, solution_fname=f"sol_self_{src_name}.h5")


    # step 9: subtract
    # img for subtract
    run_wsclean_imaging(applycal_avg_self_ms, str(output_dir / f"{output_prefix}_{src_name}_avg_self"), 
            niter=5000, mgain=0.8, horizon_mask=2, auto_pix_fov=False, size=128, scale='1arcmin', save_source_list=True, no_negative=True,
            join_channels=True, channels_out=4, fit_spectral_pol=2)

    subtract_ms = output_dir / f"{output_prefix}_{src_name}_subtract.ms"

    sourcelist_fname = output_dir / f"{output_prefix}_{src_name}_avg_self-sources.txt"
    run_dp3_subtract(applycal_self_ms, subtract_ms, sourcelist_fname)

    run_wsclean_imaging(subtract_ms, str(output_dir / f"{output_prefix}_{src_name}_subtract_self"), 
            niter=1000, mgain=0.8, horizon_mask=3, auto_pix_fov=False, size=128, scale='1arcmin', save_source_list=True, no_negative=True,
            join_channels=True, channels_out=4, fit_spectral_pol=2)

    # step 10 reset zenith with https://wsclean.readthedocs.io/en/latest/chgcentre.html
    run_reset_zenith(subtract_ms)
    run_reset_zenith(applycal_self_ms)

    # step 11: imaging
    run_wsclean_imaging(subtract_ms, str(output_dir / f"{output_prefix}_{src_name}_subtract"), 
                       niter=8000, mgain=0.8, horizon_mask=3, auto_pix_fov=False, save_source_list=False, multiscale=True, minuv_l=3)

    # step 6: imaging
    
    total_elapsed = time.time() - pipeline_start
    print("="*60)
    print(f"Peeling pipeline completed! (Total time: {total_elapsed:.1f}s)")
    print(f"Combined MS: {combined_ms}")
    print(f"Initial image: {output_dir / f'{output_prefix}_initial'}*.fits")
    print("="*60)
    
    return str(combined_ms)

def main():
    parser = argparse.ArgumentParser(
        description="LWA Peeling Pipeline: Multi-MS wideband processing",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("--ms-list", type=str, nargs='+', required=True,
                        help="List of input measurement sets")
    parser.add_argument("--bcal-list", type=str, nargs='+', required=True,
                        help="List of bandpass calibration tables (must match MS list order)")
    parser.add_argument("--output-prefix", type=str, default="peel",
                        help="Output file prefix")
    parser.add_argument("--output-dir", type=str, default=None,
                        help="Output directory (default: parent of first MS)")
    
    args = parser.parse_args()
    
    # Run the peeling pipeline
    run_peel_pipeline(args.ms_list, args.bcal_list, 
                     output_prefix=args.output_prefix,
                     output_dir=args.output_dir)

if __name__ == "__main__":
    main()