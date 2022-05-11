#!/bin/sh

#OV 1.0
python minimize.py --run_ori 410_HPK_CH09_OV1.0_Tm30_run160 --run_target 479_HPK_CH08_OV1.0_Tm40_run166 --norm 5.85 --scale 0.83 --smear 240
python minimize.py --run_ori 410_HPK_CH09_OV1.0_Tm30_run160 --run_target 479_HPK_CH10_OV1.0_Tm40_run170 --norm 5.85 --scale 0.83 --smear 240

#OV 1.2
python minimize.py --run_ori 479_HPK_CH10_OV1.0_Tm40_run170 --run_target 479_HPK_CH08_OV1.2_Tm40_run167 --norm 6.4 --scale 0.765 --smear 370
python minimize.py --run_ori 479_HPK_CH10_OV1.0_Tm40_run170 --run_target 479_HPK_CH10_OV1.2_Tm40_run171 --norm 6.4 --scale 0.765 --smear 370

#OV0.8
python minimize.py --run_ori 410_HPK_CH09_OV0.8_Tm30_run159 --run_target 479_HPK_CH08_OV0.8_Tm40_run165 --norm 6.2 --scale 0.84 --smear 185
