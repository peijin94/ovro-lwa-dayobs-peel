import casatasks
from ovrolwasolar import flagging
import argparse
from casatasks import flagdata, split


parser = argparse.ArgumentParser()
parser.add_argument("input_ms", type=str)
parser.add_argument("gaintable", type=str)

args = parser.parse_args()
input_ms = args.input_ms
gaintable = args.gaintable


flagging.flag_bad_ants(input_ms)

casatasks.applycal( vis=input_ms,gaintable=gaintable, applymode='calflag')