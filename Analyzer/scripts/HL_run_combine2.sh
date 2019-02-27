date1=$1
date2=$2
model=$3
energy=$4
which=$5

if [[ $which == "all" ]]; then
    
    # Combine
    cd Combine; cmsenv; cd -
    python scripts/run_combine2.py --ext_only --blind -n 5 --nproc=20 -d "syst_results/run_"$date1"_syst" -s -m $model WAna_nj45 WAna_nj6 TopAna
    python scripts/run_combine2.py --ext_only --blind -n 5 --nproc=20 -d "syst_results/run_"$date1"_syst" -s -m $model WAna_nj45
    python scripts/run_combine2.py --ext_only --blind -n 5 --nproc=20 -d "syst_results/run_"$date1"_syst" -s -m $model WAna_nj6
    python scripts/run_combine2.py --ext_only --blind -n 5 --nproc=20 -d "syst_results/run_"$date1"_syst" -s -m $model TopAna
    python scripts/run_combine2.py --ext_only --blind -n 5 --nproc=20 -d "syst_results/run_"$date1"_syst" -s -m $model WAna_nj45 WAna_nj6
    python scripts/run_combine2.py --ext_only --blind -n 5 --nproc=20 -d "syst_results/run_"$date1"_syst" -m $model WAna_nj45 WAna_nj6 TopAna
    python scripts/run_combine2.py --ext_only --blind -n 5 --nproc=20 -d "syst_results/run_"$date1"_syst" -m $model WAna_nj45
    python scripts/run_combine2.py --ext_only --blind -n 5 --nproc=20 -d "syst_results/run_"$date1"_syst" -m $model WAna_nj6
    python scripts/run_combine2.py --ext_only --blind -n 5 --nproc=20 -d "syst_results/run_"$date1"_syst" -m $model TopAna
    python scripts/run_combine2.py --ext_only --blind -n 5 --nproc=20 -d "syst_results/run_"$date1"_syst" -m $model WAna_nj45 WAna_nj6
    cmsenv
    
    python scripts/Get2DContour.py -e $energy -d "syst_results/run_"$date1"_syst" -m $model -s -b WAna_nj45
    python scripts/Get2DContour.py -e $energy -d "syst_results/run_"$date1"_syst" -m $model -s -b WAna_nj6
    python scripts/Get2DContour.py -e $energy -d "syst_results/run_"$date1"_syst" -m $model -s -b TopAna
    python scripts/Get2DContour.py -e $energy -d "syst_results/run_"$date1"_syst" -m $model -s -b WAna_nj45_WAna_nj6
    python scripts/Get2DContour.py -e $energy -d "syst_results/run_"$date1"_syst" -m $model -s -b WAna_nj45_WAna_nj6_TopAna
    python scripts/Get2DContour.py -e $energy -d "syst_results/run_"$date1"_syst" -m $model --blind -b WAna_nj45
    python scripts/Get2DContour.py -e $energy -d "syst_results/run_"$date1"_syst" -m $model --blind -b WAna_nj6
    python scripts/Get2DContour.py -e $energy -d "syst_results/run_"$date1"_syst" -m $model --blind -b TopAna
    python scripts/Get2DContour.py -e $energy -d "syst_results/run_"$date1"_syst" -m $model --blind -b WAna_nj45_WAna_nj6
    python scripts/Get2DContour.py -e $energy -d "syst_results/run_"$date1"_syst" -m $model --blind -b WAna_nj45_WAna_nj6_TopAna
    
elif [[ $which == "final" ]]; then
    
    # Combine
    cd Combine; cmsenv; cd -
    python scripts/run_combine2.py --ext_only --blind -n 5 --nproc=20 -d "syst_results/run_"$date1"_syst" -s -m $model WAna_nj45 WAna_nj6 TopAna
    python scripts/run_combine2.py --ext_only --blind -n 5 --nproc=20 -d "syst_results/run_"$date1"_syst" -m $model WAna_nj45 WAna_nj6 TopAna
    cmsenv
    
    python scripts/Get2DContour.py -e $energy -d "syst_results/run_"$date1"_syst" -m $model -s -b WAna_nj45_WAna_nj6_TopAna
    python scripts/Get2DContour.py -e $energy -d "syst_results/run_"$date1"_syst" -m $model --blind -b WAna_nj45_WAna_nj6_TopAna
    
elif [[ $which == "rest" ]]; then
    
    # Combine
    cd Combine; cmsenv; cd -
    #python scripts/run_combine2.py --ext_only --blind -n 5 --nproc=20 -d "syst_results/run_"$date1"_syst" -s -m $model WAna_nj45
    #python scripts/run_combine2.py --ext_only --blind -n 5 --nproc=20 -d "syst_results/run_"$date1"_syst" -s -m $model WAna_nj6
    python scripts/run_combine2.py --ext_only --blind -n 5 --nproc=20 -d "syst_results/run_"$date1"_syst" -s -m $model TopAna
    #python scripts/run_combine2.py --ext_only --blind -n 5 --nproc=20 -d "syst_results/run_"$date1"_syst" -s -m $model WAna_nj45 WAna_nj6
    #python scripts/run_combine2.py --ext_only --blind -n 5 --nproc=20 -d "syst_results/run_"$date1"_syst" -m $model WAna_nj45
    #python scripts/run_combine2.py --ext_only --blind -n 5 --nproc=20 -d "syst_results/run_"$date1"_syst" -m $model WAna_nj6
    python scripts/run_combine2.py --ext_only --blind -n 5 --nproc=20 -d "syst_results/run_"$date1"_syst" -m $model TopAna
    #python scripts/run_combine2.py --ext_only --blind -n 5 --nproc=20 -d "syst_results/run_"$date1"_syst" -m $model WAna_nj45 WAna_nj6
    cmsenv
    
    #python scripts/Get2DContour.py -e $energy -d "syst_results/run_"$date1"_syst" -m $model -s -b WAna_nj45
    #python scripts/Get2DContour.py -e $energy -d "syst_results/run_"$date1"_syst" -m $model -s -b WAna_nj6
    python scripts/Get2DContour.py -e $energy -d "syst_results/run_"$date1"_syst" -m $model -s -b TopAna
    #python scripts/Get2DContour.py -e $energy -d "syst_results/run_"$date1"_syst" -m $model -s -b WAna_nj45_WAna_nj6
    #python scripts/Get2DContour.py -e $energy -d "syst_results/run_"$date1"_syst" -m $model --blind -b WAna_nj45
    #python scripts/Get2DContour.py -e $energy -d "syst_results/run_"$date1"_syst" -m $model --blind -b WAna_nj6
    python scripts/Get2DContour.py -e $energy -d "syst_results/run_"$date1"_syst" -m $model --blind -b TopAna
    #python scripts/Get2DContour.py -e $energy -d "syst_results/run_"$date1"_syst" -m $model --blind -b WAna_nj45_WAna_nj6
    
fi
