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
import config
from source_list import get_time_mjd, get_Sun_RA_DEC, mask_far_Sun_sources

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

def run_combine_ms(input_ms_list, output_ms):
    """DP3 combine multiple MS files"""
    print(f"Step : DP3 combine MS - {len(input_ms_list)} files -> {output_ms}")
    start_time = time.time()
    
    output_path = Path(output_ms)
    
    # Create msin list string: [ms1,ms2,ms3,...]
    msin_str = "[" + ",".join([str(Path(ms)) for ms in input_ms_list]) + "]"
    
    parset_content = f"""msin={msin_str}
        msout={str(output_path)}
        msin.datacolumn=CORRECTED_DATA
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
        applied_bp_ms_list.append(str(ms))
    
    # Step 2: Combine MS files
    combined_ms = output_dir / f"{timestamp}_combined.ms"
    print(f"\nCombining {len(applied_bp_ms_list)} MS files...")
    run_combine_ms(applied_bp_ms_list, str(combined_ms))
    
    # Step 3: Initial imaging
    print(f"\nGenerating initial image...")
    run_wsclean_imaging(combined_ms, str(output_dir / f"{output_prefix}_initial"), 
                       niter=8000, mgain=0.8, horizon_mask=3, auto_pix_fov=False)
    
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