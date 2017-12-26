date1=$1
date2=$2
model=$3

python scripts/run_syst.py -d "results/run_"$date1"_syst" -m $model -b WAna_nj35
python scripts/run_syst.py -d "results/run_"$date1"_syst" -m $model -b WAna_nj45
python scripts/run_syst.py -d "results/run_"$date1"_syst" -m $model -b WAna_nj46
python scripts/run_syst.py -d "results/run_"$date1"_syst" -m $model -b WAna_nj6
python scripts/run_syst.py -d "results/run_"$date1"_syst" -m $model -b WAna_nj7
python scripts/run_syst.py -d "results/run_"$date2"_syst_TopAna" -m $model -b TopAna
cp -p "syst_results/run_"$date2"_syst_TopAna/cards/RazorBoost_SMS-"$model"_"*"_TopAna."* "syst_results/run_"$date1"_syst/cards/"
cd Combine; cmsenv; cd -
python scripts/run_combine.py -n 5 --nproc=4 -d "syst_results/run_"$date1"_syst" -m $model WAna_nj35
python scripts/run_combine.py -n 5 --nproc=4 -d "syst_results/run_"$date1"_syst" -m $model WAna_nj45
python scripts/run_combine.py -n 5 --nproc=4 -d "syst_results/run_"$date1"_syst" -m $model WAna_nj46
python scripts/run_combine.py -n 5 --nproc=4 -d "syst_results/run_"$date1"_syst" -m $model WAna_nj6
python scripts/run_combine.py -n 5 --nproc=4 -d "syst_results/run_"$date1"_syst" -m $model WAna_nj7
python scripts/run_combine.py -n 5 --nproc=4 -d "syst_results/run_"$date1"_syst" -m $model TopAna
python scripts/run_combine.py -n 5 --nproc=4 -d "syst_results/run_"$date1"_syst" -m $model WAna_nj35 WAna_nj6
python scripts/run_combine.py -n 5 --nproc=4 -d "syst_results/run_"$date1"_syst" -m $model WAna_nj45 WAna_nj6
python scripts/run_combine.py -n 5 --nproc=4 -d "syst_results/run_"$date1"_syst" -m $model WAna_nj46 WAna_nj7
python scripts/run_combine.py -n 5 --nproc=4 -d "syst_results/run_"$date1"_syst" -m $model WAna_nj35 WAna_nj6 TopAna
python scripts/run_combine.py -n 5 --nproc=4 -d "syst_results/run_"$date1"_syst" -m $model WAna_nj45 WAna_nj6 TopAna
python scripts/run_combine.py -n 5 --nproc=4 -d "syst_results/run_"$date1"_syst" -m $model WAna_nj46 WAna_nj7 TopAna
cmsenv
python scripts/Get2DContour.py -d "syst_results/run_"$date1"_syst" -m $model -b WAna_nj35
python scripts/Get2DContour.py -d "syst_results/run_"$date1"_syst" -m $model -b WAna_nj45
python scripts/Get2DContour.py -d "syst_results/run_"$date1"_syst" -m $model -b WAna_nj46
python scripts/Get2DContour.py -d "syst_results/run_"$date1"_syst" -m $model -b WAna_nj6
python scripts/Get2DContour.py -d "syst_results/run_"$date1"_syst" -m $model -b WAna_nj7
python scripts/Get2DContour.py -d "syst_results/run_"$date1"_syst" -m $model -b TopAna
python scripts/Get2DContour.py -d "syst_results/run_"$date1"_syst" -m $model -b WAna_nj35_WAna_nj6
python scripts/Get2DContour.py -d "syst_results/run_"$date1"_syst" -m $model -b WAna_nj45_WAna_nj6
python scripts/Get2DContour.py -d "syst_results/run_"$date1"_syst" -m $model -b WAna_nj46_WAna_nj7
python scripts/Get2DContour.py -d "syst_results/run_"$date1"_syst" -m $model -b WAna_nj35_WAna_nj6_TopAna
python scripts/Get2DContour.py -d "syst_results/run_"$date1"_syst" -m $model -b WAna_nj45_WAna_nj6_TopAna
python scripts/Get2DContour.py -d "syst_results/run_"$date1"_syst" -m $model -b WAna_nj46_WAna_nj7_TopAna
mkdir -p "config/"$date1
sed "s;SYST_DIR;syst_results/run_"$date1"_syst;g;s;Model;"$model";g;s;Box;WAna_nj35;g;s;boxname;Wn35;g;" config/template.cfg > "config/"$date1"/"$model"_RazorBoost_Wn35.cfg"
sed "s;SYST_DIR;syst_results/run_"$date1"_syst;g;s;Model;"$model";g;s;Box;WAna_nj45;g;s;boxname;Wn45;g;" config/template.cfg > "config/"$date1"/"$model"_RazorBoost_Wn45.cfg"
sed "s;SYST_DIR;syst_results/run_"$date1"_syst;g;s;Model;"$model";g;s;Box;WAna_nj46;g;s;boxname;Wn46;g;" config/template.cfg > "config/"$date1"/"$model"_RazorBoost_Wn46.cfg"
sed "s;SYST_DIR;syst_results/run_"$date1"_syst;g;s;Model;"$model";g;s;Box;WAna_nj6;g;s;boxname;Wn6;g;" config/template.cfg > "config/"$date1"/"$model"_RazorBoost_Wn6.cfg"
sed "s;SYST_DIR;syst_results/run_"$date1"_syst;g;s;Model;"$model";g;s;Box;WAna_nj7;g;s;boxname;Wn7;g;" config/template.cfg > "config/"$date1"/"$model"_RazorBoost_Wn7.cfg"
sed "s;SYST_DIR;syst_results/run_"$date1"_syst;g;s;Model;"$model";g;s;Box;TopAna;g;s;boxname;Top;g;" config/template.cfg > "config/"$date1"/"$model"_RazorBoost_Top.cfg"
sed "s;SYST_DIR;syst_results/run_"$date1"_syst;g;s;Model;"$model";g;s;Box;WAna_nj35_WAna_nj6;g;s;boxname;W;g;" config/template.cfg > "config/"$date1"/"$model"_RazorBoost_Wn35_Wn6.cfg"
sed "s;SYST_DIR;syst_results/run_"$date1"_syst;g;s;Model;"$model";g;s;Box;WAna_nj45_WAna_nj6;g;s;boxname;W;g;" config/template.cfg > "config/"$date1"/"$model"_RazorBoost_Wn45_Wn6.cfg"
sed "s;SYST_DIR;syst_results/run_"$date1"_syst;g;s;Model;"$model";g;s;Box;WAna_nj46_WAna_nj7;g;s;boxname;W;g;" config/template.cfg > "config/"$date1"/"$model"_RazorBoost_Wn46_Wn7.cfg"
sed "s;SYST_DIR;syst_results/run_"$date1"_syst;g;s;Model;$model;g;s;Box;WAna_nj35_WAna_nj6_TopAna;g;s;boxname;Boost;g;" config/template.cfg > "config/"$date1"/"$model"_RazorBoost_Wn35_Wn6_Top.cfg"
sed "s;SYST_DIR;syst_results/run_"$date1"_syst;g;s;Model;$model;g;s;Box;WAna_nj45_WAna_nj6_TopAna;g;s;boxname;Boost;g;" config/template.cfg > "config/"$date1"/"$model"_RazorBoost_Wn45_Wn6_Top.cfg"
sed "s;SYST_DIR;syst_results/run_"$date1"_syst;g;s;Model;$model;g;s;Box;WAna_nj46_WAna_nj7_TopAna;g;s;boxname;Boost;g;" config/template.cfg > "config/"$date1"/"$model"_RazorBoost_Wn46_Wn7_Top.cfg"
python scripts/makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Wn35.cfg" $model"_Wn35_"
python scripts/makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Wn45.cfg" $model"_Wn45_"
python scripts/makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Wn46.cfg" $model"_Wn46_"
python scripts/makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Wn6.cfg"  $model"_Wn6_"
python scripts/makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Wn7.cfg"  $model"_Wn7_"
python scripts/makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Top.cfg" $model"_Top_"
python scripts/makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Wn35_Wn6.cfg" $model"_Wn35_Wn6_"
python scripts/makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Wn45_Wn6.cfg" $model"_Wn45_Wn6_"
python scripts/makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Wn46_Wn7.cfg" $model"_Wn46_Wn7_"
python scripts/makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Wn35_Wn6_Top.cfg" $model"_Wn35_Wn6_Top_"
python scripts/makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Wn45_Wn6_Top.cfg" $model"_Wn45_Wn6_Top_"
python scripts/makeSMSplots.py "config/"$date1"/"$model"_RazorBoost_Wn46_Wn7_Top.cfg" $model"_Wn46_Wn7_Top_"
mkdir -p "Plots/limits/"$date1
mv "$model"*"_XSEC."* "Plots/limits/"$date1"/"
