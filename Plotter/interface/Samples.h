#ifndef Samples_h
#define Samples_h

#include <vector>
#include <string>

class Samples {
  
public:
  
  Samples() { NSample_=0; }
  ~Samples() {}
  
private:
  
  size_t NSample_;
  std::vector<std::string> PFNames_;
  std::vector<std::string> LatexNames_;
  std::vector<std::string> Colors_;
  std::vector<std::vector<std::string> > DirNames_;
  std::vector<std::vector<double> > ScaleFactors_;
  std::vector<size_t > DirToIndex_;
  
public:
  
  typedef struct SubSample { std::string dir; double k_factor; } SubSample;
  void AddSample(std::string name, std::string latex, std::string col, std::vector<SubSample> sample) {
    PFNames_.push_back(name);
    LatexNames_.push_back(latex);
    Colors_.push_back(col);
    DirNames_.push_back(std::vector<std::string>());
    KFactors_.push_back(std::vector<double>());
    for ( SubSample subsample : sample ) {
      DirToIndex_.push_back(NSample_);
      DirNames_[DirNames_.size()-1].push_back(subsample.dir);
      KFactors_[KFactors_.size()-1].push_back(subsample.k_factor);
    }
    ++NSample_;
  }
  
  std::vector<std::string> GetListOfDirectories() {
    std::vector<std::string> dirs;
    for ( auto subdirs : DirNames_ ) for ( auto subdir : subdirs ) dirs.push_back(subdir);
    return dirs;
  }
  
  std::vector<double> GetKFactors() {
    std::vector<double> v_k_factor;
    for ( auto subsample_k_factors : KFactors_ ) for ( auto k : subsample_k_factors) v_k_factor.push_back(k);
    return v_k_factor;
  }
  
  std::vector<size_t > GetDirToIndex() { return DirToIndex_; }
  
  std::string GetPFNames() {
    std::string names;
    for ( std::string name : PFNames_ ) names += name+";";
    return names.substr(0, names.size()-1);
  }
    
  std::string GetLatexNames() {
    std::string names;
    for ( std::string name : LatexNames_ ) names += name+";";
    return names.substr(0, names.size()-1);
  }
    
  std::string GetColors() {
    std::string colors;
    for ( std::string col : Colors_ ) colors += col+",";
    return colors.substr(0, colors.size()-1);
  }
  
  std::vector<std::string> GetListOfNames() { return PFNames_; }
  std::vector<std::string> GetListOfLat() { return LatexNames_; }
};

#endif
