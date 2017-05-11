// Example To Call Generate Designer From C++

// Include Files
#include "rt_nonfinite.h"
#include "internal_design_filter_cg.h"
#include "ex_main.h"
#include <stdio.h>
#include <string>

// Function Prototype Declarations
static void main_internal_design_filter_cg();

// Function Definitions
static void main_internal_design_filter_cg()
{

  // Initialize input arguments
  double Rdata = 7680000;
  double Fpass = 2250000;
  double Fstop = 2750000;
  double caldiv = 32;
  double FIR = 2;
  double HB1 = 2;
  double DAC_div = 1;

  std::string SType = "Lowpass";
  const char *Type = SType.c_str();
  std::string SRxTx = "Tx";
  const char *RxTx = SRxTx.c_str();

  double RFbw = 4236204;
  double converter_rate = 122880000;
  double PLL_rate = 983040000;
  double Fcenter = 0;
  double wnom = 4400000;
  double FIRdBmin = 0;
  double int_FIR = 1;
  double PLL_mult = 8;
  double Apass = 0.1250;
  double Astop = 85;
  double phEQ = -1;
  double HB2 = 2;
  double HB3 = 3;

  // Initialize output
  short outputTaps[128];

  // Call the entry-point 'internal_design_filter_cg'.
  internal_design_filter_cg(Rdata, Fpass, Fstop, caldiv, FIR, HB1, PLL_mult,
    Apass, Astop, phEQ, HB2, HB3, Type, RxTx, RFbw, DAC_div, converter_rate,
    PLL_rate, Fcenter, wnom, FIRdBmin, int_FIR, outputTaps);

  // Show taps
  for (int k=0; k<128; k++)
    printf("%i\n",outputTaps[k]);
}


int main(int, const char * const [])
{
  // Initialize the application.
  internal_design_filter_cg_initialize();

  // Invoke the entry-point functions.
  main_internal_design_filter_cg();

  // Terminate the application.
  internal_design_filter_cg_terminate();
  return 0;
}

// [EOF]
