date1=$1
date2=$2
model=$3
energy=$4
lumi=$5
scenario=$6

# Cards
#python scripts/HL_run_syst.py -s $scenario -l $lumi -e $energy -d "results/run_"$date1"_syst" -m $model -b WAna_nj45
#python scripts/HL_run_syst.py -s $scenario -l $lumi -e $energy -d "results/run_"$date1"_syst" -m $model -b WAna_nj6
python scripts/HL_run_syst.py -s $scenario -l $lumi -e $energy -d "results/run_"$date2"_syst_TopAna" -m $model -b TopAna
cp -p "syst_results/run_"$date2"_syst_TopAna/cards/RazorBoost_SMS-"$model"_"*"_TopAna."* "syst_results/run_"$date1"_syst/cards/"
