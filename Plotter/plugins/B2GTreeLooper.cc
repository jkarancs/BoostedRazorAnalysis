#ifndef B2GTreeLooper_h
#define B2GTreeLooper_h

#include "TFile.h"
#include "TH3.h"
#include "TStopwatch.h"
#include "TChain.h"
#include <string>
#include <iomanip>
#include <iostream>
#include <sstream>

class B2GTreeLooper {
public:
  B2GTreeLooper(int nthstat = 1, bool show_progress = 0) :
    nthstat_(nthstat),show_progress_(show_progress) {
    init_();
  }
  ~B2GTreeLooper() {
    delete sw_;
  }
  
private:
  const int nthstat_;
  const bool show_progress_;
  
  void init_() {
    if (show_progress_) {
      sw_ = new TStopwatch;
      processed_entries_=0;
      progress_=0;
      step_size_ = 0.0001;
      total_entries_=0;
    }
    it_sample=-1;
  }
  
  // For B2GTreeLooper
  std::vector<TChain*> samples_;
  TObjArray* files_;
  TFile* file_;
  Long64_t entry_;
  Long64_t all_entries_;
  const char* treename_;
  
  // For ProgressEstimator
  TStopwatch* sw_;
  Long64_t total_entries_;
  Long64_t processed_entries_;
  double progress_;
  double step_size_;
  

  void ProgressEstimator_() {
    if (processed_entries_==0) sw_->Start(0);
    if (double(processed_entries_)/double(total_entries_)>progress_+step_size_) {
      progress_ += step_size_;
      sw_->Stop();
      double totaltime = sw_->RealTime();
      sw_->Continue();
      double entriespersec = double(processed_entries_)/totaltime;
      double timeleft = double(total_entries_-processed_entries_)/ entriespersec;
      int time_h = int(timeleft/3600);
      int time_m = (timeleft-(time_h*3600))/60;
      int time_s = timeleft-(time_h*3600)-(time_m)*60;
      step_size_ = ((time_h)?30:1)*(entriespersec/double(total_entries_));
      std::cout<<"\r"<<std::flush;
      printf("----------- %03.1f%% ( %02d:%02d:%02d Left ) -----------", progress_*100, time_h, time_m, time_s);
    }
    processed_entries_++;
  }
  
public:
  size_t it_sample, it_file;
  std::vector<Long64_t > nevents;
  
  //Modifiers
  void AddFile(std::string fileaddress, int new_loop=1) {
    std::string treename = "/B2GTTreeMaker/B2GTree";
    size_t size = samples_.size();
    if (new_loop!=0) {
      std::stringstream ss;
      ss<<"_"<<size+1;
      samples_.push_back(new TChain((std::string("filechain")+ss.str()).c_str()));
      size = samples_.size();
    }
    if (new_loop!=-1) samples_[size-1]->Add((fileaddress+treename).c_str());

    // Calculate total number of entries - For ProgressEstimator
    Long64_t sample_entries = 0;
    if (nthstat_==1) {
      sample_entries = samples_[size-1]->GetEntries();
    } else for (int nf=0; nf<samples_[size-1]->GetListOfFiles()->GetEntries(); ++nf) {
      TFile f(samples_[size-1]->GetListOfFiles()->At(nf)->GetTitle());
      sample_entries += ((TTree*)f.Get(samples_[size-1]->GetListOfFiles()->At(nf)->GetName()))->GetEntries()/nthstat_;
    }
    total_entries_+=sample_entries;
    step_size_ = 50000/(total_entries_==0 ? 1 : total_entries_);
    nevents.push_back(sample_entries);
  }
  
  bool LoopOnSamples() {
    ++it_sample;
    if (it_sample<samples_.size()) {
      files_ = samples_[it_sample]->GetListOfFiles();
      it_file=0;
      return 1;
    } else {
      return 0;
    }
  }
  
  bool LoopOnFiles() {
    if (it_file<(size_t)files_->GetEntries()) {
      //std::cout<<"Opening file: "<<files_->At(it_file)->GetTitle()<<std::endl;
      file_ = TFile::Open(files_->At(it_file)->GetTitle());
      treename_ = files_->At(it_file)->GetName();
      all_entries_ = ((TTree*)file_->Get(treename_))->GetEntries()/nthstat_;
      if (!show_progress_) std::cout<<"Opened file: "<<files_->At(it_file)->GetTitle()<<std::endl;
      ++it_file;
      entry_=-1;
      return 1;
    } else {
      return 0;
    }
  }
  
  bool LoopOnEntries() { 
    ++entry_; 
    if (show_progress_) ProgressEstimator_();
    return entry_<all_entries_;
  }
  
  // Accessors
  const char* TreeName() { return treename_; }
  
  TFile* CurrentFile() { return file_; }
  
  Long64_t CurrentEntry() { return entry_; }


};

#endif
