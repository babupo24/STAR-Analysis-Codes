//
// File generated by rootcint at Thu May 12 11:45:00 2022

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME dOsl73_gcc485dIobjdIStRootdIStSpinPooldIStFmsEventDisplaydIStFmsEventDisplay_Cint
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "StFmsEventDisplay_Cint.h"

#include "TClass.h"
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

// Direct notice to TROOT of the dictionary's loading.
namespace {
   static struct DictInit {
      DictInit() {
         ROOT::RegisterModule();
      }
   } __TheDictionaryInitializer;
}

// START OF SHADOWS

namespace ROOTShadow {
   namespace Shadow {
   } // of namespace Shadow
} // of namespace ROOTShadow
// END OF SHADOWS

namespace ROOTDict {
   void StFmsEventDisplay_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_StFmsEventDisplay(void *p = 0);
   static void *newArray_StFmsEventDisplay(Long_t size, void *p);
   static void delete_StFmsEventDisplay(void *p);
   static void deleteArray_StFmsEventDisplay(void *p);
   static void destruct_StFmsEventDisplay(void *p);

   // Function generating the singleton type initializer
   static ROOT::TGenericClassInfo *GenerateInitInstanceLocal(const ::StFmsEventDisplay*)
   {
      ::StFmsEventDisplay *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::StFmsEventDisplay >(0);
      static ::ROOT::TGenericClassInfo 
         instance("StFmsEventDisplay", ::StFmsEventDisplay::Class_Version(), ".sl73_gcc485/obj/StRoot/StSpinPool/StFmsEventDisplay/StFmsEventDisplay.h", 30,
                  typeid(::StFmsEventDisplay), ::ROOT::DefineBehavior(ptr, ptr),
                  &::StFmsEventDisplay::Dictionary, isa_proxy, 4,
                  sizeof(::StFmsEventDisplay) );
      instance.SetNew(&new_StFmsEventDisplay);
      instance.SetNewArray(&newArray_StFmsEventDisplay);
      instance.SetDelete(&delete_StFmsEventDisplay);
      instance.SetDeleteArray(&deleteArray_StFmsEventDisplay);
      instance.SetDestructor(&destruct_StFmsEventDisplay);
      return &instance;
   }
   ROOT::TGenericClassInfo *GenerateInitInstance(const ::StFmsEventDisplay*)
   {
      return GenerateInitInstanceLocal((::StFmsEventDisplay*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::StFmsEventDisplay*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOTDict

//______________________________________________________________________________
atomic_TClass_ptr StFmsEventDisplay::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *StFmsEventDisplay::Class_Name()
{
   return "StFmsEventDisplay";
}

//______________________________________________________________________________
const char *StFmsEventDisplay::ImplFileName()
{
   return ::ROOTDict::GenerateInitInstanceLocal((const ::StFmsEventDisplay*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int StFmsEventDisplay::ImplFileLine()
{
   return ::ROOTDict::GenerateInitInstanceLocal((const ::StFmsEventDisplay*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void StFmsEventDisplay::Dictionary()
{
   fgIsA = ::ROOTDict::GenerateInitInstanceLocal((const ::StFmsEventDisplay*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *StFmsEventDisplay::Class()
{
   if (!fgIsA) { R__LOCKGUARD2(gCINTMutex); if(!fgIsA) {fgIsA = ::ROOTDict::GenerateInitInstanceLocal((const ::StFmsEventDisplay*)0x0)->GetClass();} }
   return fgIsA;
}

//______________________________________________________________________________
void StFmsEventDisplay::Streamer(TBuffer &R__b)
{
   // Stream an object of class StFmsEventDisplay.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(StFmsEventDisplay::Class(),this);
   } else {
      R__b.WriteClassBuffer(StFmsEventDisplay::Class(),this);
   }
}

//______________________________________________________________________________
void StFmsEventDisplay::ShowMembers(TMemberInspector &R__insp)
{
      // Inspect the data members of an object of class StFmsEventDisplay.
      TClass *R__cl = ::StFmsEventDisplay::IsA();
      if (R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*mFmsDbMaker", &mFmsDbMaker);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*mFmsColl", &mFmsColl);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*mApplication", &mApplication);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*mCanvas", &mCanvas);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "mFilter", &mFilter);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "mNEvents", &mNEvents);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "mNAccepted", &mNAccepted);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "mMaxEvents", &mMaxEvents);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*mFilename", &mFilename);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*mFile", &mFile);
      StMaker::ShowMembers(R__insp);
}

namespace ROOTDict {
   // Wrappers around operator new
   static void *new_StFmsEventDisplay(void *p) {
      return  p ? new(p) ::StFmsEventDisplay : new ::StFmsEventDisplay;
   }
   static void *newArray_StFmsEventDisplay(Long_t nElements, void *p) {
      return p ? new(p) ::StFmsEventDisplay[nElements] : new ::StFmsEventDisplay[nElements];
   }
   // Wrapper around operator delete
   static void delete_StFmsEventDisplay(void *p) {
      delete ((::StFmsEventDisplay*)p);
   }
   static void deleteArray_StFmsEventDisplay(void *p) {
      delete [] ((::StFmsEventDisplay*)p);
   }
   static void destruct_StFmsEventDisplay(void *p) {
      typedef ::StFmsEventDisplay current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOTDict for class ::StFmsEventDisplay

/********************************************************
* .sl73_gcc485/obj/StRoot/StSpinPool/StFmsEventDisplay/StFmsEventDisplay_Cint.cxx
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************/

#ifdef G__MEMTEST
#undef malloc
#undef free
#endif

#if defined(__GNUC__) && __GNUC__ >= 4 && ((__GNUC_MINOR__ == 2 && __GNUC_PATCHLEVEL__ >= 1) || (__GNUC_MINOR__ >= 3))
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

extern "C" void G__cpp_reset_tagtableStFmsEventDisplay_Cint();

extern "C" void G__set_cpp_environmentStFmsEventDisplay_Cint() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("StFmsEventDisplay.h");
  G__cpp_reset_tagtableStFmsEventDisplay_Cint();
}
#include <new>
extern "C" int G__cpp_dllrevStFmsEventDisplay_Cint() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* StFmsEventDisplay */
static int G__StFmsEventDisplay_Cint_617_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   StFmsEventDisplay* p = NULL;
   char* gvp = (char*) G__getgvp();
   switch (libp->paran) {
   case 1:
     //m: 1
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new StFmsEventDisplay((const Char_t*) G__int(libp->para[0]));
     } else {
       p = new((void*) gvp) StFmsEventDisplay((const Char_t*) G__int(libp->para[0]));
     }
     break;
   case 0:
     int n = G__getaryconstruct();
     if (n) {
       if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
         p = new StFmsEventDisplay[n];
       } else {
         p = new((void*) gvp) StFmsEventDisplay[n];
       }
     } else {
       if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
         p = new StFmsEventDisplay;
       } else {
         p = new((void*) gvp) StFmsEventDisplay;
       }
     }
     break;
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_StFmsEventDisplay));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__StFmsEventDisplay_Cint_617_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((StFmsEventDisplay*) G__getstructoffset())->setFilter((int) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__StFmsEventDisplay_Cint_617_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((StFmsEventDisplay*) G__getstructoffset())->setMaxEvents((int) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__StFmsEventDisplay_Cint_617_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((StFmsEventDisplay*) G__getstructoffset())->setFileName((char*) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__StFmsEventDisplay_Cint_617_0_9(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) StFmsEventDisplay::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__StFmsEventDisplay_Cint_617_0_10(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) StFmsEventDisplay::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__StFmsEventDisplay_Cint_617_0_11(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) StFmsEventDisplay::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__StFmsEventDisplay_Cint_617_0_12(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      StFmsEventDisplay::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__StFmsEventDisplay_Cint_617_0_16(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((StFmsEventDisplay*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__StFmsEventDisplay_Cint_617_0_17(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) StFmsEventDisplay::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__StFmsEventDisplay_Cint_617_0_18(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) StFmsEventDisplay::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__StFmsEventDisplay_Cint_617_0_19(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) StFmsEventDisplay::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__StFmsEventDisplay_Cint_617_0_20(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) StFmsEventDisplay::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__StFmsEventDisplay_Cint_617_0_21(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   StFmsEventDisplay* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new StFmsEventDisplay(*(StFmsEventDisplay*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_StFmsEventDisplay));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef StFmsEventDisplay G__TStFmsEventDisplay;
static int G__StFmsEventDisplay_Cint_617_0_22(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   char* gvp = (char*) G__getgvp();
   long soff = G__getstructoffset();
   int n = G__getaryconstruct();
   //
   //has_a_delete: 1
   //has_own_delete1arg: 0
   //has_own_delete2arg: 0
   //
   if (!soff) {
     return(1);
   }
   if (n) {
     if (gvp == (char*)G__PVOID) {
       delete[] (StFmsEventDisplay*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((StFmsEventDisplay*) (soff+(sizeof(StFmsEventDisplay)*i)))->~G__TStFmsEventDisplay();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (StFmsEventDisplay*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((StFmsEventDisplay*) (soff))->~G__TStFmsEventDisplay();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* StFmsEventDisplay */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncStFmsEventDisplay_Cint {
 public:
  G__Sizep2memfuncStFmsEventDisplay_Cint(): p(&G__Sizep2memfuncStFmsEventDisplay_Cint::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncStFmsEventDisplay_Cint::*p)();
};

size_t G__get_sizep2memfuncStFmsEventDisplay_Cint()
{
  G__Sizep2memfuncStFmsEventDisplay_Cint a;
  G__setsizep2memfunc((int)a.sizep2memfunc());
  return((size_t)a.sizep2memfunc());
}


/*********************************************************
* virtual base class offset calculation interface
*********************************************************/

   /* Setting up class inheritance */

/*********************************************************
* Inheritance information setup/
*********************************************************/
extern "C" void G__cpp_setup_inheritanceStFmsEventDisplay_Cint() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_StFmsEventDisplay))) {
     StFmsEventDisplay *G__Lderived;
     G__Lderived=(StFmsEventDisplay*)0x1000;
     {
       StMaker *G__Lpbase=(StMaker*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_StFmsEventDisplay),G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_StMaker),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
     {
       TDataSet *G__Lpbase=(TDataSet*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_StFmsEventDisplay),G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_TDataSet),(long)G__Lpbase-(long)G__Lderived,1,0);
     }
     {
       TNamed *G__Lpbase=(TNamed*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_StFmsEventDisplay),G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_TNamed),(long)G__Lpbase-(long)G__Lderived,1,0);
     }
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_StFmsEventDisplay),G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,0);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableStFmsEventDisplay_Cint() {

   /* Setting up typedef entry */
   G__search_typename2("Char_t",99,-1,0,-1);
   G__setnewtype(-1,"Signed Character 1 byte (char)",0);
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<std::bidirectional_iterator_tag,TObject*,std::ptrdiff_t,const TObject**,const TObject*&>",117,G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*,std::ptrdiff_t,const TObject**,const TObject*&>",117,G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*>",117,G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*,long>",117,G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*,long,const TObject**>",117,G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("map<std::string,TObjArray*>",117,G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_maplEstringcOTObjArraymUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOTObjArraymUgRsPgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("map<string,TObjArray*>",117,G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_maplEstringcOTObjArraymUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOTObjArraymUgRsPgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("map<string,TObjArray*>",117,G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_maplEstringcOTObjArraymUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOTObjArraymUgRsPgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("map<string,TObjArray*,less<string> >",117,G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_maplEstringcOTObjArraymUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOTObjArraymUgRsPgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TVectorT<Float_t>",117,G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_TVectorTlEfloatgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TVectorT<Double_t>",117,G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_TVectorTlEdoublegR),0,-1);
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* StFmsEventDisplay */
static void G__setup_memvarStFmsEventDisplay(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_StFmsEventDisplay));
   { StFmsEventDisplay *p; p=(StFmsEventDisplay*)0x1000; if (p) { }
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_StFmsDbMaker),-1,-1,4,"mFmsDbMaker=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_StFmsCollection),-1,-1,4,"mFmsColl=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_TApplication),-1,-1,4,"mApplication=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_TCanvas),-1,-1,4,"mCanvas=",0,(char*)NULL);
   G__memvar_setup((void*)0,105,0,0,-1,-1,-1,4,"mFilter=",0,(char*)NULL);
   G__memvar_setup((void*)0,105,0,0,-1,-1,-1,4,"mNEvents=",0,(char*)NULL);
   G__memvar_setup((void*)0,105,0,0,-1,-1,-1,4,"mNAccepted=",0,(char*)NULL);
   G__memvar_setup((void*)0,105,0,0,-1,-1,-1,4,"mMaxEvents=",0,(char*)NULL);
   G__memvar_setup((void*)0,67,0,0,-1,-1,-1,4,"mFilename=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_TFile),-1,-1,4,"mFile=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_TClass),G__defined_typename("atomic_TClass_ptr"),-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarStFmsEventDisplay_Cint() {
}
/***********************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
***********************************************************/

/*********************************************************
* Member function information setup for each class
*********************************************************/
static void G__setup_memfuncStFmsEventDisplay(void) {
   /* StFmsEventDisplay */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_StFmsEventDisplay));
   G__memfunc_setup("StFmsEventDisplay",1733,G__StFmsEventDisplay_Cint_617_0_1, 105, G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_StFmsEventDisplay), -1, 0, 1, 1, 1, 0, "C - 'Char_t' 10 '\"FmsED\"' name", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Init",404,(G__InterfaceMethod) NULL,105, -1, G__defined_typename("Int_t"), 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Make",382,(G__InterfaceMethod) NULL,105, -1, G__defined_typename("Int_t"), 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Finish",609,(G__InterfaceMethod) NULL,105, -1, G__defined_typename("Int_t"), 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("setFilter",946,G__StFmsEventDisplay_Cint_617_0_5, 121, -1, -1, 0, 1, 1, 1, 0, "i - - 0 - v", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("setMaxEvents",1255,G__StFmsEventDisplay_Cint_617_0_6, 121, -1, -1, 0, 1, 1, 1, 0, "i - - 0 - v", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("setFileName",1101,G__StFmsEventDisplay_Cint_617_0_7, 121, -1, -1, 0, 1, 1, 1, 0, "C - - 0 - file", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetCVS",524,(G__InterfaceMethod) NULL, 67, -1, -1, 0, 0, 1, 4, 9, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Class",502,G__StFmsEventDisplay_Cint_617_0_9, 85, G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&StFmsEventDisplay::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__StFmsEventDisplay_Cint_617_0_10, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&StFmsEventDisplay::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__StFmsEventDisplay_Cint_617_0_11, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&StFmsEventDisplay::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__StFmsEventDisplay_Cint_617_0_12, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&StFmsEventDisplay::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,(G__InterfaceMethod) NULL,85, G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TMemberInspector' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__StFmsEventDisplay_Cint_617_0_16, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - ClassDef_StreamerNVirtual_b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__StFmsEventDisplay_Cint_617_0_17, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&StFmsEventDisplay::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__StFmsEventDisplay_Cint_617_0_18, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&StFmsEventDisplay::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__StFmsEventDisplay_Cint_617_0_19, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&StFmsEventDisplay::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__StFmsEventDisplay_Cint_617_0_20, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&StFmsEventDisplay::DeclFileLine) ), 0);
   // automatic copy constructor
   G__memfunc_setup("StFmsEventDisplay", 1733, G__StFmsEventDisplay_Cint_617_0_21, (int) ('i'), G__get_linked_tagnum(&G__StFmsEventDisplay_CintLN_StFmsEventDisplay), -1, 0, 1, 1, 1, 0, "u 'StFmsEventDisplay' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~StFmsEventDisplay", 1859, G__StFmsEventDisplay_Cint_617_0_22, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 1);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncStFmsEventDisplay_Cint() {
}

/*********************************************************
* Global variable information setup for each class
*********************************************************/
static void G__cpp_setup_global0() {

   /* Setting up global variables */
   G__resetplocal();

}

static void G__cpp_setup_global1() {
}

static void G__cpp_setup_global2() {
}

static void G__cpp_setup_global3() {
}

static void G__cpp_setup_global4() {
}

static void G__cpp_setup_global5() {
}

static void G__cpp_setup_global6() {
}

static void G__cpp_setup_global7() {
}

static void G__cpp_setup_global8() {

   G__resetglobalenv();
}
extern "C" void G__cpp_setup_globalStFmsEventDisplay_Cint() {
  G__cpp_setup_global0();
  G__cpp_setup_global1();
  G__cpp_setup_global2();
  G__cpp_setup_global3();
  G__cpp_setup_global4();
  G__cpp_setup_global5();
  G__cpp_setup_global6();
  G__cpp_setup_global7();
  G__cpp_setup_global8();
}

/*********************************************************
* Global function information setup for each class
*********************************************************/
static void G__cpp_setup_func0() {
   G__lastifuncposition();

}

static void G__cpp_setup_func1() {
}

static void G__cpp_setup_func2() {
}

static void G__cpp_setup_func3() {
}

static void G__cpp_setup_func4() {
}

static void G__cpp_setup_func5() {
}

static void G__cpp_setup_func6() {
}

static void G__cpp_setup_func7() {
}

static void G__cpp_setup_func8() {
}

static void G__cpp_setup_func9() {
}

static void G__cpp_setup_func10() {
}

static void G__cpp_setup_func11() {
}

static void G__cpp_setup_func12() {
}

static void G__cpp_setup_func13() {
}

static void G__cpp_setup_func14() {
}

static void G__cpp_setup_func15() {
}

static void G__cpp_setup_func16() {
}

static void G__cpp_setup_func17() {
}

static void G__cpp_setup_func18() {
}

static void G__cpp_setup_func19() {
}

static void G__cpp_setup_func20() {
}

static void G__cpp_setup_func21() {
}

static void G__cpp_setup_func22() {
}

static void G__cpp_setup_func23() {
}

static void G__cpp_setup_func24() {
}

static void G__cpp_setup_func25() {
}

static void G__cpp_setup_func26() {
}

static void G__cpp_setup_func27() {

   G__resetifuncposition();
}

extern "C" void G__cpp_setup_funcStFmsEventDisplay_Cint() {
  G__cpp_setup_func0();
  G__cpp_setup_func1();
  G__cpp_setup_func2();
  G__cpp_setup_func3();
  G__cpp_setup_func4();
  G__cpp_setup_func5();
  G__cpp_setup_func6();
  G__cpp_setup_func7();
  G__cpp_setup_func8();
  G__cpp_setup_func9();
  G__cpp_setup_func10();
  G__cpp_setup_func11();
  G__cpp_setup_func12();
  G__cpp_setup_func13();
  G__cpp_setup_func14();
  G__cpp_setup_func15();
  G__cpp_setup_func16();
  G__cpp_setup_func17();
  G__cpp_setup_func18();
  G__cpp_setup_func19();
  G__cpp_setup_func20();
  G__cpp_setup_func21();
  G__cpp_setup_func22();
  G__cpp_setup_func23();
  G__cpp_setup_func24();
  G__cpp_setup_func25();
  G__cpp_setup_func26();
  G__cpp_setup_func27();
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__StFmsEventDisplay_CintLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__StFmsEventDisplay_CintLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__StFmsEventDisplay_CintLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__StFmsEventDisplay_CintLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__StFmsEventDisplay_CintLN_TNamed = { "TNamed" , 99 , -1 };
G__linked_taginfo G__StFmsEventDisplay_CintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__StFmsEventDisplay_CintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__StFmsEventDisplay_CintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__StFmsEventDisplay_CintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__StFmsEventDisplay_CintLN_TDataSet = { "TDataSet" , 99 , -1 };
G__linked_taginfo G__StFmsEventDisplay_CintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR = { "iterator<bidirectional_iterator_tag,TObject*,long,const TObject**,const TObject*&>" , 115 , -1 };
G__linked_taginfo G__StFmsEventDisplay_CintLN_maplEstringcOTObjArraymUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOTObjArraymUgRsPgRsPgR = { "map<string,TObjArray*,less<string>,allocator<pair<const string,TObjArray*> > >" , 99 , -1 };
G__linked_taginfo G__StFmsEventDisplay_CintLN_TVectorTlEfloatgR = { "TVectorT<float>" , 99 , -1 };
G__linked_taginfo G__StFmsEventDisplay_CintLN_TVectorTlEdoublegR = { "TVectorT<double>" , 99 , -1 };
G__linked_taginfo G__StFmsEventDisplay_CintLN_TFile = { "TFile" , 99 , -1 };
G__linked_taginfo G__StFmsEventDisplay_CintLN_StMaker = { "StMaker" , 99 , -1 };
G__linked_taginfo G__StFmsEventDisplay_CintLN_TCanvas = { "TCanvas" , 99 , -1 };
G__linked_taginfo G__StFmsEventDisplay_CintLN_TApplication = { "TApplication" , 99 , -1 };
G__linked_taginfo G__StFmsEventDisplay_CintLN_StFmsDbMaker = { "StFmsDbMaker" , 99 , -1 };
G__linked_taginfo G__StFmsEventDisplay_CintLN_StFmsCollection = { "StFmsCollection" , 99 , -1 };
G__linked_taginfo G__StFmsEventDisplay_CintLN_StFmsEventDisplay = { "StFmsEventDisplay" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableStFmsEventDisplay_Cint() {
  G__StFmsEventDisplay_CintLN_TClass.tagnum = -1 ;
  G__StFmsEventDisplay_CintLN_TBuffer.tagnum = -1 ;
  G__StFmsEventDisplay_CintLN_TMemberInspector.tagnum = -1 ;
  G__StFmsEventDisplay_CintLN_TObject.tagnum = -1 ;
  G__StFmsEventDisplay_CintLN_TNamed.tagnum = -1 ;
  G__StFmsEventDisplay_CintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__StFmsEventDisplay_CintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__StFmsEventDisplay_CintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__StFmsEventDisplay_CintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__StFmsEventDisplay_CintLN_TDataSet.tagnum = -1 ;
  G__StFmsEventDisplay_CintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR.tagnum = -1 ;
  G__StFmsEventDisplay_CintLN_maplEstringcOTObjArraymUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOTObjArraymUgRsPgRsPgR.tagnum = -1 ;
  G__StFmsEventDisplay_CintLN_TVectorTlEfloatgR.tagnum = -1 ;
  G__StFmsEventDisplay_CintLN_TVectorTlEdoublegR.tagnum = -1 ;
  G__StFmsEventDisplay_CintLN_TFile.tagnum = -1 ;
  G__StFmsEventDisplay_CintLN_StMaker.tagnum = -1 ;
  G__StFmsEventDisplay_CintLN_TCanvas.tagnum = -1 ;
  G__StFmsEventDisplay_CintLN_TApplication.tagnum = -1 ;
  G__StFmsEventDisplay_CintLN_StFmsDbMaker.tagnum = -1 ;
  G__StFmsEventDisplay_CintLN_StFmsCollection.tagnum = -1 ;
  G__StFmsEventDisplay_CintLN_StFmsEventDisplay.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableStFmsEventDisplay_Cint() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__StFmsEventDisplay_CintLN_TClass);
   G__get_linked_tagnum_fwd(&G__StFmsEventDisplay_CintLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__StFmsEventDisplay_CintLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__StFmsEventDisplay_CintLN_TObject);
   G__get_linked_tagnum_fwd(&G__StFmsEventDisplay_CintLN_TNamed);
   G__get_linked_tagnum_fwd(&G__StFmsEventDisplay_CintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__StFmsEventDisplay_CintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__StFmsEventDisplay_CintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__StFmsEventDisplay_CintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__StFmsEventDisplay_CintLN_TDataSet);
   G__get_linked_tagnum_fwd(&G__StFmsEventDisplay_CintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR);
   G__get_linked_tagnum_fwd(&G__StFmsEventDisplay_CintLN_maplEstringcOTObjArraymUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOTObjArraymUgRsPgRsPgR);
   G__get_linked_tagnum_fwd(&G__StFmsEventDisplay_CintLN_TVectorTlEfloatgR);
   G__get_linked_tagnum_fwd(&G__StFmsEventDisplay_CintLN_TVectorTlEdoublegR);
   G__get_linked_tagnum_fwd(&G__StFmsEventDisplay_CintLN_TFile);
   G__get_linked_tagnum_fwd(&G__StFmsEventDisplay_CintLN_StMaker);
   G__get_linked_tagnum_fwd(&G__StFmsEventDisplay_CintLN_TCanvas);
   G__get_linked_tagnum_fwd(&G__StFmsEventDisplay_CintLN_TApplication);
   G__get_linked_tagnum_fwd(&G__StFmsEventDisplay_CintLN_StFmsDbMaker);
   G__get_linked_tagnum_fwd(&G__StFmsEventDisplay_CintLN_StFmsCollection);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__StFmsEventDisplay_CintLN_StFmsEventDisplay),sizeof(StFmsEventDisplay),-1,324864,(char*)NULL,G__setup_memvarStFmsEventDisplay,G__setup_memfuncStFmsEventDisplay);
}
extern "C" void G__cpp_setupStFmsEventDisplay_Cint(void) {
  G__check_setup_version(30051515,"G__cpp_setupStFmsEventDisplay_Cint()");
  G__set_cpp_environmentStFmsEventDisplay_Cint();
  G__cpp_setup_tagtableStFmsEventDisplay_Cint();

  G__cpp_setup_inheritanceStFmsEventDisplay_Cint();

  G__cpp_setup_typetableStFmsEventDisplay_Cint();

  G__cpp_setup_memvarStFmsEventDisplay_Cint();

  G__cpp_setup_memfuncStFmsEventDisplay_Cint();
  G__cpp_setup_globalStFmsEventDisplay_Cint();
  G__cpp_setup_funcStFmsEventDisplay_Cint();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncStFmsEventDisplay_Cint();
  return;
}
class G__cpp_setup_initStFmsEventDisplay_Cint {
  public:
    G__cpp_setup_initStFmsEventDisplay_Cint() { G__add_setup_func("StFmsEventDisplay_Cint",(G__incsetup)(&G__cpp_setupStFmsEventDisplay_Cint)); G__call_setup_funcs(); }
   ~G__cpp_setup_initStFmsEventDisplay_Cint() { G__remove_setup_func("StFmsEventDisplay_Cint"); }
};
G__cpp_setup_initStFmsEventDisplay_Cint G__cpp_setup_initializerStFmsEventDisplay_Cint;

