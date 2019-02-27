date1=$1
date2=$2
model=$3

# Expected limits
mkdir -p "config/"$date1
#sed "s;SYST_DIR;syst_results/run_"$date1";g;s;Model;"$model";g;s;Box;WAna_nj45;g;s;boxname;W_4-5_jet;g;" config/HE2_template_exp.cfg > "config/"$date1"/"$model"_RazorBoost_Wn45_exp.cfg"
#sed "s;SYST_DIR;syst_results/run_"$date1";g;s;Model;"$model";g;s;Box;WAna_nj6;g;s;boxname;W_6#leq_jet;g;" config/HE2_template_exp.cfg > "config/"$date1"/"$model"_RazorBoost_Wn6_exp.cfg"
#sed "s;SYST_DIR;syst_results/run_"$date1";g;s;Model;"$model";g;s;Box;TopAna;g;s;boxname;Top;g;" config/HE2_template_exp.cfg > "config/"$date1"/"$model"_RazorBoost_Top_exp.cfg"
#sed "s;SYST_DIR;syst_results/run_"$date1";g;s;Model;"$model";g;s;Box;WAna_nj45_WAna_nj6;g;s;boxname;W_4-5_and_6#leq_jet;g;" config/HE2_template_exp.cfg > "config/"$date1"/"$model"_RazorBoost_Wn45_Wn6_exp.cfg"
sed "s;SYST_DIR;syst_results/run_"$date1";g;s;Model;$model;g;s;Box;WAna_nj45_WAna_nj6_TopAna;g;s;boxname;;g;" config/HE2_template_exp.cfg > "config/"$date1"/"$model"_RazorBoost_Wn45_Wn6_Top_exp.cfg"
#python scripts/HE2_makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Wn45_exp.cfg" $model"_Wn45_"
#python scripts/HE2_makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Wn6_exp.cfg"  $model"_Wn6_"
#python scripts/HE2_makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Top_exp.cfg" $model"_Top_"
#python scripts/HE2_makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Wn45_Wn6_exp.cfg" $model"_Wn45_Wn6_"
python scripts/HE2_makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Wn45_Wn6_Top_exp.cfg" $model"_Wn45_Wn6_Top_"
mkdir -p "Plots/exp_limits_try/"$date1
mv "$model"*"_XSEC."* "Plots/exp_limits_try/"$date1"/"

# Significance
mkdir -p "config/"$date1
#sed "s;SYST_DIR;syst_results/run_"$date1";g;s;Model;"$model";g;s;Box;WAna_nj45;g;s;boxname;W_4-5_jet;g;" config/HE2_template_sig.cfg > "config/"$date1"/"$model"_RazorBoost_Wn45_sig.cfg"
#sed "s;SYST_DIR;syst_results/run_"$date1";g;s;Model;"$model";g;s;Box;WAna_nj6;g;s;boxname;W_6#leq_jet;g;" config/HE2_template_sig.cfg > "config/"$date1"/"$model"_RazorBoost_Wn6_sig.cfg"
#sed "s;SYST_DIR;syst_results/run_"$date1";g;s;Model;"$model";g;s;Box;TopAna;g;s;boxname;Top;g;" config/HE2_template_sig.cfg > "config/"$date1"/"$model"_RazorBoost_Top_sig.cfg"
#sed "s;SYST_DIR;syst_results/run_"$date1";g;s;Model;"$model";g;s;Box;WAna_nj45_WAna_nj6;g;s;boxname;W_4-5_and_6#leq_jet;g;" config/HE2_template_sig.cfg > "config/"$date1"/"$model"_RazorBoost_Wn45_Wn6_sig.cfg"
sed "s;SYST_DIR;syst_results/run_"$date1";g;s;Model;$model;g;s;Box;WAna_nj45_WAna_nj6_TopAna;g;s;boxname;;g;" config/HE2_template_sig.cfg > "config/"$date1"/"$model"_RazorBoost_Wn45_Wn6_Top_sig.cfg"
#python scripts/HE2_makeSMSplots_sig.py "config/"$date1"/"$model"_RazorBoost_Wn45_sig.cfg" $model"_Wn45_"
#python scripts/HE2_makeSMSplots_sig.py "config/"$date1"/"$model"_RazorBoost_Wn6_sig.cfg"  $model"_Wn6_"
#python scripts/HE2_makeSMSplots_sig.py "config/"$date1"/"$model"_RazorBoost_Top_sig.cfg" $model"_Top_"
#python scripts/HE2_makeSMSplots_sig.py "config/"$date1"/"$model"_RazorBoost_Wn45_Wn6_sig.cfg" $model"_Wn45_Wn6_"
python scripts/HE2_makeSMSplots_sig.py "config/"$date1"/"$model"_RazorBoost_Wn45_Wn6_Top_sig.cfg" $model"_Wn45_Wn6_Top_"
mkdir -p "Plots/significance_try/"$date1
mv "$model"*"_SIG."* "Plots/significance_try/"$date1"/"
