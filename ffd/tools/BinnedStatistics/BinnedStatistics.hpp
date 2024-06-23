#ifndef FFD_BINNED_STAT_PACKAGE_NOT_BEEN
#define FFD_BINNED_STAT_PACKAGE_NOT_BEEN


#include"../../core.hpp"
#include<functional>
#include<map>


template<typename Field>
struct BinnedStatistics{
  long NumberOfStoredValues = 0;
  long BinWidth = 1;
  int NumberDiscardedBins=0;
  vector<Field> Data;
  vector<Real> Data2;
  Real MaximalFluctuation = 2;

  
  BinnedStatistics(){}
  

  BinnedStatistics(int MaximalNumberOfBins, int NumberDiscardedBins0=0, int BinWidth0 = 1): NumberDiscardedBins(NumberDiscardedBins0), BinWidth(BinWidth0) {
    Data = vector<Field>(MaximalNumberOfBins, 0.);
    Data2 = vector<Real>(MaximalNumberOfBins, 0.);
  }
  
  
  void StoreValue(Field value){
    int which_bin = NumberOfStoredValues/BinWidth;
    Data[which_bin] += value;
    Data2[which_bin] += pow(abs(value), 2);
    NumberOfStoredValues++;
    if(NumberOfStoredValues == BinWidth*Data.size()){
      for(int u=0; u<Data.size()/2; ++u){
	Data[u] = Data[2*u]+Data[2*u+1];
	Data2[u] = Data2[2*u]+Data2[2*u+1];
      }
      for(int u=Data.size()/2; u<Data.size(); ++u){
	Data[u] = 0;
	Data2[u] = 0;
      }
      BinWidth *= 2;
    }
  }

  
  Field Average(){
    Field ret = 0;
    long NumberOfUsableValues = NumberOfStoredValues - NumberDiscardedBins*BinWidth;
    if(NumberOfUsableValues <= 0){
      return 0;
    }
    for(int u=NumberDiscardedBins; u<Data.size(); ++u){
      ret += Data[u];
    }
    return ret/NumberOfUsableValues;
  }

  
  
  Real AbsoluteError(){
    Field avg = 0;
    Real avg2 = 0;
    int NumberOfFullBins = NumberOfStoredValues/BinWidth-NumberDiscardedBins;
    if(NumberOfFullBins <= 1){
      return MaximalFluctuation;
    }
    for(int u=NumberDiscardedBins; u<NumberOfFullBins+NumberDiscardedBins; ++u){
      avg += Data[u]/BinWidth;
      avg2 += pow(abs(Data[u]/BinWidth), 2);
    }
    avg /= NumberOfFullBins;
    avg2 /= NumberOfFullBins;
    avg2 -= pow(abs(avg), 2);
    avg2 = sqrt(avg2/(NumberOfFullBins-1));
    return avg2 + MaximalFluctuation/NumberOfStoredValues;
  }
  

  
  Real RelativeError(){
    return AbsoluteError()/Average();
  }

  
  Real AutoCorrelationTime(){
    Field avg = 0;
    Real avg2 = 0;
    long NumberOfUsableValues = NumberOfStoredValues - NumberDiscardedBins*BinWidth;
    for(int u=NumberDiscardedBins; u<Data.size(); ++u){
      avg += Data[u];
      avg2 += Data2[u];
    }
    avg /= NumberOfUsableValues;
    avg2 /= NumberOfUsableValues;
    avg2 -= pow(abs(avg), 2);
    avg2 = avg2/(NumberOfUsableValues-1);
    return .5*( pow(Error(), 2) / avg2 - 1 );
  }
};




#endif
