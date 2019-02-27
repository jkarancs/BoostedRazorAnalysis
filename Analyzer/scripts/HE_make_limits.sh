date1=$1
date2=$2
model=$3
intlumi=$4

# Expected limits
mkdir -p "config/"$date1
sed "s;SYST_DIR;syst_results/run_"$date1";g;s;Model;"$model";g;s;Box;WAna_nj45;g;s;boxname;W_4-5_jet;g;s;INTLUMI;"$intlumi"000000;g;" config/HE_template_exp.cfg > "config/"$date1"/"$model"_RazorBoost_Wn45_exp.cfg"
sed "s;SYST_DIR;syst_results/run_"$date1";g;s;Model;"$model";g;s;Box;WAna_nj6;g;s;boxname;W_6#leq_jet;g;s;INTLUMI;"$intlumi"000000;g;" config/HE_template_exp.cfg > "config/"$date1"/"$model"_RazorBoost_Wn6_exp.cfg"
sed "s;SYST_DIR;syst_results/run_"$date1";g;s;Model;"$model";g;s;Box;TopAna;g;s;boxname;Top;g;s;INTLUMI;"$intlumi"000000;g;" config/HE_template_exp.cfg > "config/"$date1"/"$model"_RazorBoost_Top_exp.cfg"
sed "s;SYST_DIR;syst_results/run_"$date1";g;s;Model;"$model";g;s;Box;WAna_nj45_WAna_nj6;g;s;boxname;W_4-5/6#leq_jet;g;s;INTLUMI;"$intlumi"000000;g;" config/HE_template_exp.cfg > "config/"$date1"/"$model"_RazorBoost_Wn45_Wn6_exp.cfg"
sed "s;SYST_DIR;syst_results/run_"$date1";g;s;Model;$model;g;s;Box;WAna_nj45_WAna_nj6_TopAna;g;s;boxname;;g;s;INTLUMI;"$intlumi"000000;g;" config/HE_template_exp.cfg > "config/"$date1"/"$model"_RazorBoost_Wn45_Wn6_Top_exp.cfg"
python scripts/HE_makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Wn45_exp.cfg" $model"_Wn45_"
python scripts/HE_makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Wn6_exp.cfg"  $model"_Wn6_"
python scripts/HE_makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Top_exp.cfg" $model"_Top_"
python scripts/HE_makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Wn45_Wn6_exp.cfg" $model"_Wn45_Wn6_"
python scripts/HE_makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Wn45_Wn6_Top_exp.cfg" $model"_Wn45_Wn6_Top_"
mkdir -p "Plots/exp_limits/"$date1
mv "$model"*"_XSEC."* "Plots/exp_limits/"$date1

# Expected limits (statonly)
mkdir -p "config/"$date1
sed "s;SYST_DIR;syst_results/run_"$date1";g;s;Model;"$model";g;s;Box;WAna_nj45;g;s;boxname;W_4-5_jet;g;s;INTLUMI;"$intlumi"000000;g;" config/HE_template_exp_statonly.cfg > "config/"$date1"/"$model"_RazorBoost_Wn45_exp_statonly.cfg"
sed "s;SYST_DIR;syst_results/run_"$date1";g;s;Model;"$model";g;s;Box;WAna_nj6;g;s;boxname;W_6#leq_jet;g;s;INTLUMI;"$intlumi"000000;g;" config/HE_template_exp_statonly.cfg > "config/"$date1"/"$model"_RazorBoost_Wn6_exp_statonly.cfg"
sed "s;SYST_DIR;syst_results/run_"$date1";g;s;Model;"$model";g;s;Box;TopAna;g;s;boxname;Top;g;s;INTLUMI;"$intlumi"000000;g;" config/HE_template_exp_statonly.cfg > "config/"$date1"/"$model"_RazorBoost_Top_exp_statonly.cfg"
sed "s;SYST_DIR;syst_results/run_"$date1";g;s;Model;"$model";g;s;Box;WAna_nj45_WAna_nj6;g;s;boxname;W_4-5/6#leq_jet;g;s;INTLUMI;"$intlumi"000000;g;" config/HE_template_exp_statonly.cfg > "config/"$date1"/"$model"_RazorBoost_Wn45_Wn6_exp_statonly.cfg"
sed "s;SYST_DIR;syst_results/run_"$date1";g;s;Model;$model;g;s;Box;WAna_nj45_WAna_nj6_TopAna;g;s;boxname;;g;s;INTLUMI;"$intlumi"000000;g;" config/HE_template_exp_statonly.cfg > "config/"$date1"/"$model"_RazorBoost_Wn45_Wn6_Top_exp_statonly.cfg"
python scripts/HE_makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Wn45_exp_statonly.cfg" $model"_Wn45_"
python scripts/HE_makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Wn6_exp_statonly.cfg"  $model"_Wn6_"
python scripts/HE_makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Top_exp_statonly.cfg" $model"_Top_"
python scripts/HE_makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Wn45_Wn6_exp_statonly.cfg" $model"_Wn45_Wn6_"
python scripts/HE_makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Wn45_Wn6_Top_exp_statonly.cfg" $model"_Wn45_Wn6_Top_"
mkdir -p "Plots/exp_limits/"$date1"_statonly"
mv "$model"*"_XSEC."* "Plots/exp_limits/"$date1"_statonly/"

# Significance (logz - manually set: sms.py, smsPlot3SIG.py)
mkdir -p "config/"$date1
sed "s;SYST_DIR;syst_results/run_"$date1";g;s;Model;"$model";g;s;Box;WAna_nj45;g;s;boxname;W_4-5_jet;g;s;INTLUMI;"$intlumi"000000;g;" config/HE_template_sig.cfg > "config/"$date1"/"$model"_RazorBoost_Wn45_sig.cfg"
sed "s;SYST_DIR;syst_results/run_"$date1";g;s;Model;"$model";g;s;Box;WAna_nj6;g;s;boxname;W_6#leq_jet;g;s;INTLUMI;"$intlumi"000000;g;" config/HE_template_sig.cfg > "config/"$date1"/"$model"_RazorBoost_Wn6_sig.cfg"
sed "s;SYST_DIR;syst_results/run_"$date1";g;s;Model;"$model";g;s;Box;TopAna;g;s;boxname;Top;g;s;INTLUMI;"$intlumi"000000;g;" config/HE_template_sig.cfg > "config/"$date1"/"$model"_RazorBoost_Top_sig.cfg"
sed "s;SYST_DIR;syst_results/run_"$date1";g;s;Model;"$model";g;s;Box;WAna_nj45_WAna_nj6;g;s;boxname;W_4-5/6#leq_jet;g;s;INTLUMI;"$intlumi"000000;g;" config/HE_template_sig.cfg > "config/"$date1"/"$model"_RazorBoost_Wn45_Wn6_sig.cfg"
sed "s;SYST_DIR;syst_results/run_"$date1";g;s;Model;$model;g;s;Box;WAna_nj45_WAna_nj6_TopAna;g;s;boxname;;g;s;INTLUMI;"$intlumi"000000;g;" config/HE_template_sig.cfg > "config/"$date1"/"$model"_RazorBoost_Wn45_Wn6_Top_sig.cfg"
python scripts/HE_makeSMSplots_sig.py "config/"$date1"/"$model"_RazorBoost_Wn45_sig.cfg" $model"_Wn45_logz_"
python scripts/HE_makeSMSplots_sig.py "config/"$date1"/"$model"_RazorBoost_Wn6_sig.cfg"  $model"_Wn6_logz_"
python scripts/HE_makeSMSplots_sig.py "config/"$date1"/"$model"_RazorBoost_Top_sig.cfg" $model"_Top_logz_"
python scripts/HE_makeSMSplots_sig.py "config/"$date1"/"$model"_RazorBoost_Wn45_Wn6_sig.cfg" $model"_Wn45_Wn6_logz_"
python scripts/HE_makeSMSplots_sig.py "config/"$date1"/"$model"_RazorBoost_Wn45_Wn6_Top_sig.cfg" $model"_Wn45_Wn6_Top_logz_"
mkdir -p "Plots/significance/"$date1
mv "$model"*"_SIG."* "Plots/significance/"$date1"/"
