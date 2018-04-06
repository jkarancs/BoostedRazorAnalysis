date1=$1
date2=$2
model=$3

# Cards
python scripts/run_syst.py -d "results/run_"$date1"_syst" -m $model -b WAna_nj45
python scripts/run_syst.py -d "results/run_"$date1"_syst" -m $model -b WAna_nj6
python scripts/run_syst.py -d "results/run_"$date2"_syst_TopAna" -m $model -b TopAna
cp -p "syst_results/run_"$date2"_syst_TopAna/cards/RazorBoost_SMS-"$model"_"*"_TopAna."* "syst_results/run_"$date1"_syst/cards/"

# Combine
cd Combine; cmsenv; cd -
python scripts/run_combine.py -n 5 --nproc=10 -d "syst_results/run_"$date1"_syst" -m $model WAna_nj45 WAna_nj6 TopAna
python scripts/run_combine.py -n 5 --nproc=10 -d "syst_results/run_"$date1"_syst" -m $model WAna_nj45
python scripts/run_combine.py -n 5 --nproc=10 -d "syst_results/run_"$date1"_syst" -m $model WAna_nj6
python scripts/run_combine.py -n 5 --nproc=10 -d "syst_results/run_"$date1"_syst" -m $model TopAna
python scripts/run_combine.py -n 5 --nproc=10 -d "syst_results/run_"$date1"_syst" -m $model WAna_nj45 WAna_nj6
cmsenv
python scripts/Get2DContour.py -d "syst_results/run_"$date1"_syst" -m $model -b WAna_nj45
python scripts/Get2DContour.py -d "syst_results/run_"$date1"_syst" -m $model -b WAna_nj6
python scripts/Get2DContour.py -d "syst_results/run_"$date1"_syst" -m $model -b TopAna
python scripts/Get2DContour.py -d "syst_results/run_"$date1"_syst" -m $model -b WAna_nj45_WAna_nj6
python scripts/Get2DContour.py -d "syst_results/run_"$date1"_syst" -m $model -b WAna_nj45_WAna_nj6_TopAna

# Expected limits
mkdir -p "config/"$date1
sed "s;SYST_DIR;syst_results/run_"$date1"_syst;g;s;Model;"$model";g;s;Box;WAna_nj45;g;s;boxname;Wn45;g;" config/template_exp.cfg > "config/"$date1"/"$model"_RazorBoost_Wn45_exp.cfg"
sed "s;SYST_DIR;syst_results/run_"$date1"_syst;g;s;Model;"$model";g;s;Box;WAna_nj6;g;s;boxname;Wn6;g;" config/template_exp.cfg > "config/"$date1"/"$model"_RazorBoost_Wn6_exp.cfg"
sed "s;SYST_DIR;syst_results/run_"$date1"_syst;g;s;Model;"$model";g;s;Box;TopAna;g;s;boxname;Top;g;" config/template_exp.cfg > "config/"$date1"/"$model"_RazorBoost_Top_exp.cfg"
sed "s;SYST_DIR;syst_results/run_"$date1"_syst;g;s;Model;"$model";g;s;Box;WAna_nj45_WAna_nj6;g;s;boxname;W;g;" config/template_exp.cfg > "config/"$date1"/"$model"_RazorBoost_Wn45_Wn6_exp.cfg"
sed "s;SYST_DIR;syst_results/run_"$date1"_syst;g;s;Model;$model;g;s;Box;WAna_nj45_WAna_nj6_TopAna;g;s;boxname;Boost;g;" config/template_exp.cfg > "config/"$date1"/"$model"_RazorBoost_Wn45_Wn6_Top_exp.cfg"
python scripts/makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Wn45_exp.cfg" $model"_Wn45_"
python scripts/makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Wn6_exp.cfg"  $model"_Wn6_"
python scripts/makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Top_exp.cfg" $model"_Top_"
python scripts/makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Wn45_Wn6_exp.cfg" $model"_Wn45_Wn6_"
python scripts/makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Wn45_Wn6_Top_exp.cfg" $model"_Wn45_Wn6_Top_"
mkdir -p "Plots/exp_limits/"$date1
mv "$model"*"_XSEC."* "Plots/exp_limits/"$date1"/"

# Observed limits
mkdir -p "config/"$date1
sed "s;SYST_DIR;syst_results/run_"$date1"_syst;g;s;Model;"$model";g;s;Box;WAna_nj45;g;s;boxname;Wn45;g;" config/template_obs.cfg > "config/"$date1"/"$model"_RazorBoost_Wn45_obs.cfg"
sed "s;SYST_DIR;syst_results/run_"$date1"_syst;g;s;Model;"$model";g;s;Box;WAna_nj6;g;s;boxname;Wn6;g;" config/template_obs.cfg > "config/"$date1"/"$model"_RazorBoost_Wn6_obs.cfg"
sed "s;SYST_DIR;syst_results/run_"$date1"_syst;g;s;Model;"$model";g;s;Box;TopAna;g;s;boxname;Top;g;" config/template_obs.cfg > "config/"$date1"/"$model"_RazorBoost_Top_obs.cfg"
sed "s;SYST_DIR;syst_results/run_"$date1"_syst;g;s;Model;"$model";g;s;Box;WAna_nj45_WAna_nj6;g;s;boxname;W;g;" config/template_obs.cfg > "config/"$date1"/"$model"_RazorBoost_Wn45_Wn6_obs.cfg"
sed "s;SYST_DIR;syst_results/run_"$date1"_syst;g;s;Model;$model;g;s;Box;WAna_nj45_WAna_nj6_TopAna;g;s;boxname;Boost;g;" config/template_obs.cfg > "config/"$date1"/"$model"_RazorBoost_Wn45_Wn6_Top_obs.cfg"
python scripts/makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Wn45_obs.cfg" $model"_Wn45_"
python scripts/makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Wn6_obs.cfg"  $model"_Wn6_"
python scripts/makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Top_obs.cfg" $model"_Top_"
python scripts/makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Wn45_Wn6_obs.cfg" $model"_Wn45_Wn6_"
python scripts/makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Wn45_Wn6_Top_obs.cfg" $model"_Wn45_Wn6_Top_"
mkdir -p "Plots/obs_limits/"$date1
mv "$model"*"_XSEC."* "Plots/obs_limits/"$date1"/"
