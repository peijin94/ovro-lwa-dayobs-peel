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
from source_list import get_time_mjd, get_Sun_RA_DEC, mask_far_Sun_sources

PIPELINE_SCRIPT_DIR = Path(__file__).parent
EXECUTABLE_DIR = Path(__file__).parent / "exe"

def run_casa_applycal(input_ms, output_ms, gaintable):
    """Apply CASA bandpass calibration"""
    print(f"Step : CASA applycal - {input_ms}")
    start_time = time.time()
    try:
        subprocess.run(["python3", str(EXECUTABLE_DIR / "flagant_applybp.py"), str(input_ms), str(output_ms), str(gaintable)], check=True)
        elapsed = time.time() - start_time
        print(f"✓ CASA applycal completed ({elapsed:.1f}s)")
    except subprocess.CalledProcessError as e:
        elapsed = time.time() - start_time
        print(f"✗ CASA applycal failed after {elapsed:.1f}s: {e.stdout}")
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
        avg.freqstep=4
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
        print(f"✗ WSClean imaging failed after {elapsed:.1f}s: {e}")
        sys.exit(1)

def run_gaincal(input_ms, solution_fname="solution.h5", cal_type="diagonalphase"):
    """DP3 gain calibration"""
    print(f"Step : DP3 gaincal - {input_ms}")
    start_time = time.time()
    
    input_path = Path(input_ms)
#msout = /data/{input_path.name}_cal.ms
    
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
        print(f"✗ DP3 gaincal failed after {elapsed:.1f}s: {e}")
        sys.exit(1)
    
import h5py
import numpy as np

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
    


def run_dp3_subtract(input_ms, output_ms, source_list):
    """DP3 subtract"""
    print(f"Step : DP3 subtract - {input_ms} -> {output_ms}")
    start_time = time.time()
    
    input_path = Path(input_ms)
    output_path = Path(output_ms)
    source_path = Path(source_list)
    
    parset_content = f"""msin={str(input_path)}
        showprogress=False
        verbosity="quiet"
        msout={str(output_path)}
        steps=[predict]
        predict.type=predict
        predict.sourcedb={str(source_path)}
        predict.operation=subtract
        """
        
    cmd = ["DP3", *parset_content.split()]
    
    try:
        subprocess.run(cmd, check=True)
        elapsed = time.time() - start_time
        print(f"✓ DP3 subtract completed ({elapsed:.1f}s): {output_ms}")
    except subprocess.CalledProcessError as e:
        elapsed = time.time() - start_time
        print(f"✗ DP3 subtract failed after {elapsed:.1f}s: {e}")
        sys.exit(1)
    

def phaseshift_to_sun(ms_file, output_ms):
    """Phase shift MS to Sun's coordinates using DP3 PhaseShift step."""
    ms_path = Path(ms_file)
    output_path = Path(output_ms)
    if not ms_path.exists():
        raise FileNotFoundError(f"MS file not found: {ms_path}")
    start_time = time.time()

    # Get Sun position
    time_mjd = get_time_mjd(str(ms_path))
    sun_ra, sun_dec = get_Sun_RA_DEC(time_mjd)
    
    # Create parset content - use simple filenames
    parset_content = f"""msin={str(ms_path)}
        msout={str(output_path)}
        showprogress=False
        verbosity="quiet"
        steps=[phaseshift]
        phaseshift.type=phaseshift
        phaseshift.phasecenter=[{sun_ra}deg,{sun_dec}deg]
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

def run_calib_pipeline(raw_ms, gaintable, output_prefix="proc", plot_mid_steps=False, rm_ms_tmp=False, DEBUG=False, fch_img=True, mfs_img=False):
    """Run complete processing pipeline"""
    
    pipeline_start = time.time()
    raw_path = Path(raw_ms)
    data_dir = raw_path.parent
    
    # Define intermediate file paths
    applied_bp_ms = data_dir / f"{raw_path.stem}_applied_bp.ms"
    flagged_avg_ms = data_dir / f"{raw_path.stem}_flagged_avg.ms"
    solution_file = data_dir / f"{output_prefix}_solution.h5"
    final_ms = data_dir / f"{raw_path.stem}_{output_prefix}_final.ms"
    
    print("="*60)
    print("LWA Quick Processing Pipeline")
    print("="*60)
    print(f"Input: {raw_ms}")
    print(f"Gaintable: {gaintable}")
    print(f"Output prefix: {output_prefix}")
    print("="*60)

    # Step 1: casatools applycal
    run_casa_applycal(raw_ms, str(applied_bp_ms), gaintable)
    
    # Step 2: DP3 flagging and averaging, assuming corrected data column exists
    run_dp3_flag_avg(applied_bp_ms, flagged_avg_ms, strategy_file=PIPELINE_SCRIPT_DIR / "lua" / "LWA_sun_PZ.lua")

    if rm_ms_tmp:
        shutil.rmtree(raw_ms)

    current_ms = flagged_avg_ms
    # selfcal:
    run_wsclean_imaging(current_ms, str(data_dir / f"{output_prefix}_image"), niter=800, mgain=0.9,horizon_mask=5,
        save_source_list=False, auto_mask=False, auto_threshold=False)
    run_gaincal(current_ms, solution_fname=str(solution_file), cal_type="diagonalphase")
    run_applycal_dp3(current_ms,final_ms, solution_fname=str(solution_file), cal_entry_lst=["phase"])

    if rm_ms_tmp:
        shutil.rmtree(current_ms)

    # Step 6: wsclean for source subtraction
    run_wsclean_imaging(final_ms, str(data_dir / f"{output_prefix}_image_source"), niter=1500, mgain=0.9,horizon_mask=0.1 )#, multiscale=True)
    
    # Step 7: mask far Sun sources
    time_mjd = get_time_mjd(str(final_ms))
    sun_ra, sun_dec = get_Sun_RA_DEC(time_mjd)
    mask_far_Sun_sources( data_dir / f"{output_prefix}_image_source-sources.txt" , 
        data_dir / f"{output_prefix}_image_source_masked-sources.txt", 
        sun_ra, sun_dec, distance_deg=6.0)

    # Step 8: DP3 subtract sources
    subtracted_ms = data_dir / f"{output_prefix}_image_source_masked_subtracted.ms"
    print(f"Subtracting sources from {final_ms} to {subtracted_ms}", str(data_dir / f"{output_prefix}_image_source_masked-sources.txt"))
    run_dp3_subtract(final_ms, subtracted_ms, str(data_dir / f"{output_prefix}_image_source_masked-sources.txt"))

    # step 9: phaseshift to sun
    shifted_ms = data_dir / f"{output_prefix}_image_source_sun_shifted.ms"
    print(f"Phaseshifting to sun from {subtracted_ms} to {shifted_ms}")
    phaseshift_to_sun(subtracted_ms, shifted_ms)
    if rm_ms_tmp:
        shutil.rmtree(final_ms)
        shutil.rmtree(subtracted_ms)

    # final image
#    run_wsclean_imaging(shifted_ms, str(data_dir / f"{output_prefix}_image_source_sun_shifted"), auto_pix_fov=False, 
#        niter=3000, mgain=0.8, size=512, scale='1.5arcmin', save_source_list=False, weight='briggs -0.5')
    shifted_ms_avg = data_dir / f"{output_prefix}_image_source_sun_shifted_avg.ms"
    run_dp3_avg(shifted_ms, shifted_ms_avg, freq_step=4)
    
    # make a copy
    # shifted_ms_avg_copy = data_dir / f"{output_prefix}_image_source_sun_shifted_avg_copy.ms"
    # shutil.copytree(shifted_ms_avg, shifted_ms_avg_copy)

    total_elapsed = time.time() - pipeline_start

    print("="*60)
    print(f"Pipeline completed successfully! (Total time: {total_elapsed:.1f}s)")
    print("="*60)

    default_wscleancmd = "wsclean -j 8 -mem 6 -quiet -no-dirty -no-update-model-required \
        -horizon-mask 5deg -size 512 512 -scale 1.5arcmin -weight briggs -0.5 -minuv-l 10 \
        -auto-threshold 3  -niter 6000 -mgain 0.9 -beam-fitting-size 2 -pol I "
    import shlex

    if fch_img:
        time_start = time.time()
        wscleancmd = default_wscleancmd + " -join-channels -channels-out 12 -name " + f"{output_prefix}_fch " + str(shifted_ms_avg)
        subprocess.run(shlex.split(wscleancmd), check=True, capture_output=True, text=True)
        total_elapsed = time.time() - time_start
        print(f"✓ WSClean imaging completed ({total_elapsed:.1f}s): {output_prefix}_fch*.fits")

    if mfs_img:
        time_start = time.time()
        wscleancmd = default_wscleancmd + " -name " + f"{output_prefix}_mfs " + str(shifted_ms_avg)
        subprocess.run(shlex.split(wscleancmd), check=True, capture_output=True, text=True)
        total_elapsed = time.time() - time_start
        print(f"✓ WSClean imaging completed ({total_elapsed:.1f}s): {output_prefix}_mfs*.fits")

    if DEBUG:
        run_wsclean_imaging(subtracted_ms, str(data_dir / f"{output_prefix}_image_source_masked_subtracted"), niter=5000, mgain=0.9,horizon_mask=0.1)

    if plot_mid_steps:
        from script.plot_fits import plot_fits
        plot_fits(data_dir / f"{output_prefix}_image-image.fits")
        plot_fits(data_dir / f"{output_prefix}_image_source-image.fits")
        plot_fits(data_dir / f"{output_prefix}_image_source_sun_shifted-image.fits")
        if DEBUG:
            plot_fits(data_dir / f"{output_prefix}_image_source_masked_subtracted-image.fits")
        from plot_solar_image import plot_solar_image
        plot_solar_image(data_dir / f"{output_prefix}_image_source_sun_shifted-image.fits")



def main():
    parser = argparse.ArgumentParser(
        description="LWA Quick Processing Pipeline: raw MS -> CASA applycal -> DP3 flag/avg -> wsclean -> gaincal -> applycal",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Positional arguments
    parser.add_argument("raw_ms", type=str, 
                        help="Input raw measurement set path")
    parser.add_argument("gaintable", type=str, 
                        help="Calibration table path (bandpass calibration)")
    parser.add_argument("output_prefix", type=str, nargs='?', default=None,
                        help="Output file prefix (default: derived from input MS filename)")
    
    # Optional arguments
    parser.add_argument("--keep-ms-tmp", action="store_true", default=False,
                        help="Keep temporary measurement sets (default: remove them)")
    parser.add_argument("--fch-img", action="store_true", default=False,
                        help="Generate per-channel images")
    parser.add_argument("--mfs-img", action="store_true", default=False,
                        help="Generate multi-frequency synthesis image")
    
    args = parser.parse_args()
    
    # Validate inputs
    if not Path(args.raw_ms).exists():
        print(f"Error: Raw MS not found: {args.raw_ms}")
        sys.exit(1)
    
    if not Path(args.gaintable).exists():
        print(f"Error: Gaintable not found: {args.gaintable}")
        sys.exit(1)
    
    # Determine output prefix
    if args.output_prefix is None:
        args.output_prefix = Path(args.raw_ms).stem.split('.')[0]
    
    # Run the pipeline
    run_calib_pipeline( args.raw_ms,  args.gaintable,  args.output_prefix, 
        plot_mid_steps=False,  rm_ms_tmp=not args.keep_ms_tmp,  DEBUG=False,  
        fch_img=args.fch_img,  mfs_img=args.mfs_img)

if __name__ == "__main__":
    main()