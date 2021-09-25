#!/usr/local/bin/python3
import numpy as np
from config import *
from O_q_reduced import *
from test import *
from phase_diagram import *
from penetration_depth import *
from negK import *
from Emin_numerical import *
def main():
    print("(施工中) plotting publication-quality figures =3=")
    LineWidth, FontSize, LabelSize = 1,9,8

    #config_demo_plot(LineWidth, FontSize, LabelSize)
    #rippling_config_demo(LineWidth, FontSize, LabelSize)
    #phase_diagram_plot(LineWidth, FontSize, LabelSize)/
    #meron_config_demo(LineWidth, FontSize, LabelSize)
    O_q_reduced(LineWidth, FontSize, LabelSize)
    #pdepth_plot(LineWidth, FontSize, LabelSize)
    #negK_lam_Cn_plot(LineWidth, FontSize, LabelSize)
    #negK_lam_q_plot(LineWidth, FontSize, LabelSize)
    #shape_Enneper_plot(LineWidth, FontSize, LabelSize)
    #Emin_fill_between_plot_m1(LineWidth, FontSize, LabelSize)

    print("(施工中) just plotting figures for slides")
    #director_update_demo(LineWidth, FontSize, LabelSize)
    #config3d_demo_plot(LineWidth, FontSize, LabelSize)


if __name__ == "__main__":
    main()