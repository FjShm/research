#ifndef _MOTHER_ERR_INCLUDED_
#define _MOTHER_ERR_INCLUDED_
#include <iostream>
#include <string>
using namespace std ;

#ifdef WIN32
#pragma warning (disable : 4786) 
#endif

class MotherErr		// エラー処理基本クラス
{
private:
  string msg ;
  string class_name ;
  string method_name ;
public:
  MotherErr(	const string& _msg = "",const string& _method_name = "",
		const string& _class_name = "")
    :msg(_msg),class_name(_class_name),method_name(_method_name)
  {
    if (msg == "") msg = "unknown message" ;
    if (class_name == "") class_name = "unknown class name" ;
    if (method_name == "") method_name = "unknown method name" ;
  }

  void Set(const string& _msg, 
	   const string& _method_name,
	   const string& _class_name)
  { msg = _msg; class_name = _class_name; method_name = _method_name; }

  void Msg(const string& _msg){ msg = _msg; }
  void Class(const string& _class_name){ class_name = _class_name; }
  void Method(const string& _method_name){ method_name = _method_name; }

  string Msg(void){ return msg; }
  string Class(void){ return class_name; }
  string Method(void){ return method_name; }
  string Get(void) const { 
    return (class_name + "::" + method_name + " : " + msg);
  }

};

#endif //_MOTHER_ERR_INCLUDED_


#ifndef _IOPARAMETERSET_INCLUDED_
#define _IOPARAMETERSET_INCLUDED_

#include<stdlib.h>
#include<stdio.h>
#include<fstream>
#include<string>
#include<vector>
#include<algorithm>
//#include"MotherErr.h"

#ifndef STDFUNCTIONS_ARE_NOT_IN_NAMESPACE_STD
using namespace std;
#endif

//---------------------------------------------------------------------------
class SuperParameter
{
public:
  
  class ScheduleOfAParameter{
    
  public:
    
    class StepOrTimeAndValue{
    public:
      double fromTheStepOrTheTime;
      vector<string> value;
      StepOrTimeAndValue( void ){}
      StepOrTimeAndValue( double step_or_time_, const string& val_ );
      StepOrTimeAndValue( double step_or_time_, const vector<string>& val_ );
      StepOrTimeAndValue( const StepOrTimeAndValue& tv );
      ~StepOrTimeAndValue(){}
    };
    
    ScheduleOfAParameter( void );
    ScheduleOfAParameter( const string& step_or_time_,
			  double step_or_time_val,
			  const string& val_of_param );
    ScheduleOfAParameter( const string& step_or_time_,
			  double step_or_time_val,
			  const vector<string>& val_of_param );
    ScheduleOfAParameter( const ScheduleOfAParameter& sp );
    ~ScheduleOfAParameter(void){}
    friend ostream& operator<<( ostream& os, const ScheduleOfAParameter& sp );
    
    //member
    string stepOrTime;    /*! The Only Two "Step" or "Time" are permitted.*/
    vector<StepOrTimeAndValue> timeScheduledValue;
    
  }; //class ScheduleOfAParameter
  
protected:
  
  string name;
  vector<string> val;
  double min;
  double max;
  ScheduleOfAParameter ParamSchedule;
  
public:
  
  SuperParameter(const string& name_,    const string& val_,
		 const double  min_=0.0, const double max_=0.0 )
    :name(name_),min(min_),max(max_)
  { val.clear() ; val.push_back(val_) ; }
  SuperParameter(const string& name_,    const vector<string>& val_,
		 const double  min_=0.0, const double max_=0.0 )
    :name(name_),val(val_),min(min_),max(max_){}
  ~SuperParameter(){}
  SuperParameter( const SuperParameter& sp )
    :name(sp.name),val(sp.val),min(sp.min),max(sp.max){}
  
  class SuperParameterErr : public MotherErr
  {
  public:
    SuperParameterErr(const string& _msg = "",const string& _method_name = "" )
      :MotherErr( _msg, _method_name, "SuperParameter" ){}		
  };
  
  const string& Name(void) const { return(name); }
  int           Size(void) { return val.size() ; }
  double        Value(int id=0) 
  { 
    if ( (val.size() <= id) || (id < 0) )
      throw(SuperParameterErr("Parameter Array Size Error.",
			      "Value(int id=0)")) ;
    
    return( atof( val[id].c_str() ) ); 
  }
  string& StringValue(int id=0) 
  { 
    if ( (val.size() <= id) || (id < 0) )
      throw(SuperParameterErr("Parameter Array Size Error.",
			      "StringValue(int id=0)")) ;
    
    return( val[id] ); 
  }
  double& UpperLimit(void) { return( min ); }
  double& LowerLimit(void) { return( max ); }

  //  ---- check that a given string matches the name --- 
  bool MyNameIs( const string &name_ ) const ;
  
  double& operator []( const string& s );
  const ScheduleOfAParameter& ParameterSchedule(void){ return(ParamSchedule); }
  bool SetAScheduledValue( const  string& step_or_time_,
			   double step_or_time_val,
			   const  string& value_ );
  bool SetAScheduledValue( const  string& step_or_time_,
			   double step_or_time_val,
			   const  vector<string>& value_ );
  
};// class SuperParameter 
//---------------------------------------------------------------------------
//===========================================================================
//
//                       class IOParameterSet
//
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
typedef vector<SuperParameter*>      ParameterSet;
typedef ParameterSet::const_iterator ParameterSetIterator;
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
class IOParameterSet
{
  
private:
  
  ParameterSet parameters;
  
public:
  
  class IOParameterSetErr; 
  IOParameterSet(void){}
  ~IOParameterSet(void){}
  
  IOParameterSet( const string& NameOfParameterFile );
  
  // Convert all $variables to real string values.
  void ChangeAllStringsContainsDolVariableToRealVariableStrings();
  
  // Calclate unsolved string like $(A)$(/)$(B), and return solved string.
  string Change_String_Contains_DolExpr_To_SValue( const string& S );
  
  // Convert $variables in strings to real string values.
  void ChangeStringsContainsDolVariableToRealVariableStrings( vector<string>& 
															  strings );
  
  IOParameterSet( const IOParameterSet& ec );
  
  void RegisterParameter( const string& name_, const string& val_=string("0"), 
			  double min_=0.0,     double max_=1.0    );
  void RegisterParameter( const string& name_, const vector<string>& val_, 
			  double min_=0.0,     double max_=1.0    );
  
  bool ChangeValuesAccodingToSchedule( const long unsigned int& present_step, 
				       const double&            present_time );
  
  void ChangeParameter(const string& name_,const string& s,const string& val_);
  void ChangeParameter(const string& name_,const string& s,
					   const vector<string>& val_);
  
  bool SetSchedule( const string& name,      
					const string& step_or_time, 
					double step_or_time_val, 
					const string& val );
  bool SetSchedule( const string& name,      
					const string& step_or_time, 
					double step_or_time_val, 
					const vector<string>& val );
  
  void info( const string& out_dest = "") const; 
  
  // Return double parameter value, specified by "value" or "min" or "max".
  double value( const string& name_, const string& s, int id = 0 ) const;
  
  // Return double parameter value
  double value( const string& name_, int id = 0 ) const ;
  
  // Return int parameter value
  int int_value( const string& name_, int id = 0 ) const ;
  
  // Return string parameter value
  const string& string_value( const string& name_, int id = 0 ) const ;
  
  // Return bool parameter value
  bool bool_value( const string& name_, int id = 0 ) const ;
  
  // Return Parameter Object, specified by name.
  SuperParameter& GetRegisteredParameter( const string& name_ );
  
  // Check specified parameter exists or not.
  bool Exist( const string& name_ );
  
  // return number of strings given for a parameter
  int Size( const string& name_ );
  
  friend ostream& operator << ( ostream& _os, IOParameterSet& _ec );
  
  // SuperParameter& operator []( const string& name_ )
  // { return( GetRegisteredParameter( name_ ) ); }
  
};// class IOParameterSet
//----------------------------------------------------------------------------
class IOParameterSet::IOParameterSetErr:public MotherErr
{
public:
  IOParameterSetErr( const string& _msg = "", 
		     const string& _method_name = "" )
    :MotherErr( _msg, _method_name, "IOParameterSet" ){}		
};
//----------------------------------------------------------------------------
bool Check_ExistenceOfIOPfile( const int& argc );
//----------------------------------------------------------------------------
#endif //_IOPARAMETERSET_INCLUDED_
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

