import casatasks
from ovrolwasolar import flagging
import argparse
from casatasks import flagdata, split


parser = argparse.ArgumentParser()
parser.add_argument("input_ms", type=str)
parser.add_argument("output_ms", type=str)
parser.add_argument("gaintable", type=str)

args = parser.parse_args()
input_ms = args.input_ms
output_ms = args.output_ms
gaintable = args.gaintable


flagging.flag_bad_ants(input_ms)
#flagdata(vis=input_ms, mode='manual', uvrange='0.1~10lambda', flagbackup=False)
flagdata(vis=input_ms, mode='unflag', correlation='auto')
split(vis=input_ms, outputvis=output_ms, datacolumn='data', keepflags=False)

casatasks.applycal( vis=output_ms,gaintable=gaintable, applymode='calflag')