// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PackageLibCintDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "PostAnalysis/BaseFunctions.h"
#include "PostAnalysis/Constants.h"
#include "PostAnalysis/Plotting.h"
#include "PostAnalysis/UnitTest.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *BaseFunctions_Dictionary();
   static void BaseFunctions_TClassManip(TClass*);
   static void *new_BaseFunctions(void *p = 0);
   static void *newArray_BaseFunctions(Long_t size, void *p);
   static void delete_BaseFunctions(void *p);
   static void deleteArray_BaseFunctions(void *p);
   static void destruct_BaseFunctions(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::BaseFunctions*)
   {
      ::BaseFunctions *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::BaseFunctions));
      static ::ROOT::TGenericClassInfo 
         instance("BaseFunctions", "PostAnalysis/BaseFunctions.h", 10,
                  typeid(::BaseFunctions), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &BaseFunctions_Dictionary, isa_proxy, 4,
                  sizeof(::BaseFunctions) );
      instance.SetNew(&new_BaseFunctions);
      instance.SetNewArray(&newArray_BaseFunctions);
      instance.SetDelete(&delete_BaseFunctions);
      instance.SetDeleteArray(&deleteArray_BaseFunctions);
      instance.SetDestructor(&destruct_BaseFunctions);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::BaseFunctions*)
   {
      return GenerateInitInstanceLocal((::BaseFunctions*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::BaseFunctions*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *BaseFunctions_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::BaseFunctions*)0x0)->GetClass();
      BaseFunctions_TClassManip(theClass);
   return theClass;
   }

   static void BaseFunctions_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_BaseFunctions(void *p) {
      return  p ? new(p) ::BaseFunctions : new ::BaseFunctions;
   }
   static void *newArray_BaseFunctions(Long_t nElements, void *p) {
      return p ? new(p) ::BaseFunctions[nElements] : new ::BaseFunctions[nElements];
   }
   // Wrapper around operator delete
   static void delete_BaseFunctions(void *p) {
      delete ((::BaseFunctions*)p);
   }
   static void deleteArray_BaseFunctions(void *p) {
      delete [] ((::BaseFunctions*)p);
   }
   static void destruct_BaseFunctions(void *p) {
      typedef ::BaseFunctions current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::BaseFunctions

namespace {
  void TriggerDictionaryInitialization_libPackageLib_Impl() {
    static const char* headers[] = {
"PostAnalysis/BaseFunctions.h",
"PostAnalysis/Constants.h",
"PostAnalysis/Plotting.h",
"PostAnalysis/UnitTest.h",
0
    };
    static const char* includePaths[] = {
"/home/tnom6927/CTIDE/QualificationTask/PostAnalysis/PostAnalysis",
"/home/tnom6927/CTIDE/QualificationTask/PostAnalysis/PostAnalysis",
"/home/tnom6927/CTIDE/QualificationTask/PostAnalysis/PostAnalysis",
"/home/tnom6927/CTIDE/QualificationTask/PostAnalysis/PostAnalysis",
"/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AnalysisBaseExternals/21.2.47/InstallArea/x86_64-slc6-gcc62-opt/include",
"/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AnalysisBase/21.2.47/InstallArea/x86_64-slc6-gcc62-opt/RootCore/include",
"/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AnalysisBase/21.2.47/InstallArea/x86_64-slc6-gcc62-opt/RootCore/include",
"/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AnalysisBase/21.2.47/InstallArea/x86_64-slc6-gcc62-opt/RootCore/include",
"/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AnalysisBase/21.2.47/InstallArea/x86_64-slc6-gcc62-opt/RootCore/include",
"/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AnalysisBaseExternals/21.2.47/InstallArea/x86_64-slc6-gcc62-opt/include",
"/home/tnom6927/CTIDE/QualificationTask/build/PostAnalysis/CMakeFiles/makePackageLibCintDict.SsvAE1/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libPackageLib dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$PostAnalysis/BaseFunctions.h")))  BaseFunctions;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libPackageLib dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif
#ifndef HAVE_PRETTY_FUNCTION
  #define HAVE_PRETTY_FUNCTION 1
#endif
#ifndef HAVE_64_BITS
  #define HAVE_64_BITS 1
#endif
#ifndef __IDENTIFIER_64BIT__
  #define __IDENTIFIER_64BIT__ 1
#endif
#ifndef ATLAS
  #define ATLAS 1
#endif
#ifndef ROOTCORE
  #define ROOTCORE 1
#endif
#ifndef XAOD_STANDALONE
  #define XAOD_STANDALONE 1
#endif
#ifndef XAOD_ANALYSIS
  #define XAOD_ANALYSIS 1
#endif
#ifndef ROOTCORE_RELEASE_SERIES
  #define ROOTCORE_RELEASE_SERIES 25
#endif
#ifndef PACKAGE_VERSION
  #define PACKAGE_VERSION "PostAnalysis-00-00-00"
#endif
#ifndef PACKAGE_VERSION_UQ
  #define PACKAGE_VERSION_UQ PostAnalysis-00-00-00
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "PostAnalysis/BaseFunctions.h"
#include "PostAnalysis/Constants.h"
#include "PostAnalysis/Plotting.h"
#include "PostAnalysis/UnitTest.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"BaseFunctions", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libPackageLib",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libPackageLib_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libPackageLib_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libPackageLib() {
  TriggerDictionaryInitialization_libPackageLib_Impl();
}
