/*
 * problem3.c
 *
 * Code generation for model "problem3".
 *
 * Model version              : 1.198
 * Simulink Coder version : 8.6 (R2014a) 27-Dec-2013
 * C source code generated on : Tue Feb 13 12:48:06 2018
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */
#include "problem3.h"
#include "problem3_private.h"
#include "problem3_dt.h"

/* Block signals (auto storage) */
B_problem3_T problem3_B;

/* Continuous states */
X_problem3_T problem3_X;

/* Block states (auto storage) */
DW_problem3_T problem3_DW;

/* Real-time model */
RT_MODEL_problem3_T problem3_M_;
RT_MODEL_problem3_T *const problem3_M = &problem3_M_;

/*
 * Writes out MAT-file header.  Returns success or failure.
 * Returns:
 *      0 - success
 *      1 - failure
 */
int_T rt_WriteMat4FileHeader(FILE *fp, int32_T m, int32_T n, const char *name)
{
  typedef enum { ELITTLE_ENDIAN, EBIG_ENDIAN } ByteOrder;

  int16_T one = 1;
  ByteOrder byteOrder = (*((int8_T *)&one)==1) ? ELITTLE_ENDIAN : EBIG_ENDIAN;
  int32_T type = (byteOrder == ELITTLE_ENDIAN) ? 0: 1000;
  int32_T imagf = 0;
  int32_T name_len = (int32_T)strlen(name) + 1;
  if ((fwrite(&type, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&m, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&n, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&imagf, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&name_len, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(name, sizeof(char), name_len, fp) == 0)) {
    return(1);
  } else {
    return(0);
  }
}                                      /* end rt_WriteMat4FileHeader */

/*
 * This function updates continuous states using the ODE1 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE1_IntgData *id = (ODE1_IntgData *)rtsiGetSolverData(si);
  real_T *f0 = id->f[0];
  int_T i;
  int_T nXc = 4;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);
  rtsiSetdX(si, f0);
  problem3_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; i++) {
    *x += h * f0[i];
    x++;
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void problem3_output(void)
{
  /* local block i/o variables */
  real_T rtb_FromWorkspace3[4];
  real_T rtb_Backgain;
  real_T rtb_HILReadEncoderTimebase_o1;
  real_T rtb_HILReadEncoderTimebase_o2;
  real_T rtb_HILReadEncoderTimebase_o3;
  real_T *lastU;
  real_T rtb_Gain1_idx_4;
  real_T rtb_Gain1_idx_5;
  if (rtmIsMajorTimeStep(problem3_M)) {
    /* set solver stop time */
    if (!(problem3_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&problem3_M->solverInfo,
                            ((problem3_M->Timing.clockTickH0 + 1) *
        problem3_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&problem3_M->solverInfo,
                            ((problem3_M->Timing.clockTick0 + 1) *
        problem3_M->Timing.stepSize0 + problem3_M->Timing.clockTickH0 *
        problem3_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(problem3_M)) {
    problem3_M->Timing.t[0] = rtsiGetT(&problem3_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(problem3_M)) {
    /* S-Function (hil_read_encoder_timebase_block): '<S7>/HIL Read Encoder Timebase' */

    /* S-Function Block: problem3/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_read_encoder(problem3_DW.HILReadEncoderTimebase_Task, 1,
        &problem3_DW.HILReadEncoderTimebase_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(problem3_M, _rt_error_message);
      } else {
        rtb_HILReadEncoderTimebase_o1 =
          problem3_DW.HILReadEncoderTimebase_Buffer[0];
        rtb_HILReadEncoderTimebase_o2 =
          problem3_DW.HILReadEncoderTimebase_Buffer[1];
        rtb_HILReadEncoderTimebase_o3 =
          problem3_DW.HILReadEncoderTimebase_Buffer[2];
      }
    }

    /* Gain: '<S7>/Travel: Count to rad' */
    problem3_B.TravelCounttorad = problem3_P.TravelCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o1;

    /* Gain: '<S15>/Gain' */
    problem3_B.Gain = problem3_P.Gain_Gain * problem3_B.TravelCounttorad;

    /* Sum: '<Root>/Sum3' incorporates:
     *  Constant: '<Root>/travel_offset [deg]'
     */
    problem3_B.Sum3 = problem3_P.travel_offsetdeg_Value + problem3_B.Gain;

    /* Gain: '<S3>/Gain1' */
    problem3_B.Gain1 = problem3_P.Gain1_Gain * problem3_B.Sum3;

    /* Gain: '<S7>/Pitch: Count to rad' */
    problem3_B.PitchCounttorad = problem3_P.PitchCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o2;

    /* Gain: '<S12>/Gain' */
    problem3_B.Gain_i = problem3_P.Gain_Gain_a * problem3_B.PitchCounttorad;

    /* Gain: '<S4>/Gain1' */
    problem3_B.Gain1_f = problem3_P.Gain1_Gain_m * problem3_B.Gain_i;

    /* ToFile: '<Root>/To File' */
    {
      if (!(++problem3_DW.ToFile_IWORK.Decimation % 1) &&
          (problem3_DW.ToFile_IWORK.Count*2)+1 < 100000000 ) {
        FILE *fp = (FILE *) problem3_DW.ToFile_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[2];
          problem3_DW.ToFile_IWORK.Decimation = 0;
          u[0] = problem3_M->Timing.t[1];
          u[1] = problem3_B.Gain1_f;
          if (fwrite(u, sizeof(real_T), 2, fp) != 2) {
            rtmSetErrorStatus(problem3_M, "Error writing to MAT-file pitch.mat");
            return;
          }

          if (((++problem3_DW.ToFile_IWORK.Count)*2)+1 >= 100000000) {
            (void)fprintf(stdout,
                          "*** The ToFile block will stop logging data before\n"
                          "    the simulation has ended, because it has reached\n"
                          "    the maximum number of elements (100000000)\n"
                          "    allowed in MAT-file pitch.mat.\n");
          }
        }
      }
    }

    /* ToFile: '<Root>/To File1' */
    {
      if (!(++problem3_DW.ToFile1_IWORK.Decimation % 1) &&
          (problem3_DW.ToFile1_IWORK.Count*2)+1 < 100000000 ) {
        FILE *fp = (FILE *) problem3_DW.ToFile1_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[2];
          problem3_DW.ToFile1_IWORK.Decimation = 0;
          u[0] = problem3_M->Timing.t[1];
          u[1] = problem3_B.Gain1;
          if (fwrite(u, sizeof(real_T), 2, fp) != 2) {
            rtmSetErrorStatus(problem3_M, "Error writing to MAT-file travel.mat");
            return;
          }

          if (((++problem3_DW.ToFile1_IWORK.Count)*2)+1 >= 100000000) {
            (void)fprintf(stdout,
                          "*** The ToFile block will stop logging data before\n"
                          "    the simulation has ended, because it has reached\n"
                          "    the maximum number of elements (100000000)\n"
                          "    allowed in MAT-file travel.mat.\n");
          }
        }
      }
    }
  }

  /* Gain: '<S16>/Gain' incorporates:
   *  TransferFcn: '<S7>/Travel: Transfer Fcn'
   */
  problem3_B.Gain_d = (problem3_P.TravelTransferFcn_C *
                       problem3_X.TravelTransferFcn_CSTATE +
                       problem3_P.TravelTransferFcn_D *
                       problem3_B.TravelCounttorad) * problem3_P.Gain_Gain_l;

  /* Gain: '<S13>/Gain' incorporates:
   *  TransferFcn: '<S7>/Pitch: Transfer Fcn'
   */
  problem3_B.Gain_b = (problem3_P.PitchTransferFcn_C *
                       problem3_X.PitchTransferFcn_CSTATE +
                       problem3_P.PitchTransferFcn_D *
                       problem3_B.PitchCounttorad) * problem3_P.Gain_Gain_ae;

  /* FromWorkspace: '<S8>/From Workspace3' */
  {
    real_T *pDataValues = (real_T *) problem3_DW.FromWorkspace3_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *) problem3_DW.FromWorkspace3_PWORK.TimePtr;
    int_T currTimeIndex = problem3_DW.FromWorkspace3_IWORK.PrevIndex;
    real_T t = problem3_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[140]) {
      currTimeIndex = 139;
    } else {
      if (t < pTimeValues[currTimeIndex]) {
        while (t < pTimeValues[currTimeIndex]) {
          currTimeIndex--;
        }
      } else {
        while (t >= pTimeValues[currTimeIndex + 1]) {
          currTimeIndex++;
        }
      }
    }

    problem3_DW.FromWorkspace3_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 4; ++elIdx) {
              (&rtb_FromWorkspace3[0])[elIdx] = pDataValues[currTimeIndex];
              pDataValues += 141;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 4; ++elIdx) {
              (&rtb_FromWorkspace3[0])[elIdx] = pDataValues[currTimeIndex + 1];
              pDataValues += 141;
            }
          }
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;

        {
          int_T elIdx;
          for (elIdx = 0; elIdx < 4; ++elIdx) {
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            (&rtb_FromWorkspace3[0])[elIdx] = (real_T) rtInterpolate(d1, d2, f1,
              f2);
            pDataValues += 141;
          }
        }
      }
    }
  }

  /* FromWorkspace: '<S8>/From Workspace4' */
  {
    real_T *pDataValues = (real_T *) problem3_DW.FromWorkspace4_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *) problem3_DW.FromWorkspace4_PWORK.TimePtr;
    int_T currTimeIndex = problem3_DW.FromWorkspace4_IWORK.PrevIndex;
    real_T t = problem3_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[140]) {
      currTimeIndex = 139;
    } else {
      if (t < pTimeValues[currTimeIndex]) {
        while (t < pTimeValues[currTimeIndex]) {
          currTimeIndex--;
        }
      } else {
        while (t >= pTimeValues[currTimeIndex + 1]) {
          currTimeIndex++;
        }
      }
    }

    problem3_DW.FromWorkspace4_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          rtb_Backgain = pDataValues[currTimeIndex];
        } else {
          rtb_Backgain = pDataValues[currTimeIndex + 1];
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        rtb_Backgain = (real_T) rtInterpolate(d1, d2, f1, f2);
        pDataValues += 141;
      }
    }
  }

  /* Sum: '<S8>/Sum4' incorporates:
   *  Gain: '<S5>/Gain1'
   *  Gain: '<S8>/Gain'
   *  Sum: '<S8>/Sum3'
   */
  problem3_B.Sum4 = rtb_Backgain - ((((problem3_P.Gain1_Gain_b * problem3_B.Sum3
    - rtb_FromWorkspace3[0]) * problem3_P.K[0] + (problem3_P.Gain1_Gain_b *
    problem3_B.Gain_d - rtb_FromWorkspace3[1]) * problem3_P.K[1]) +
    (problem3_P.Gain1_Gain_b * problem3_B.Gain_i - rtb_FromWorkspace3[2]) *
    problem3_P.K[2]) + (problem3_P.Gain1_Gain_b * problem3_B.Gain_b -
                        rtb_FromWorkspace3[3]) * problem3_P.K[3]);
  if (rtmIsMajorTimeStep(problem3_M)) {
    /* ToFile: '<Root>/To File2' */
    {
      if (!(++problem3_DW.ToFile2_IWORK.Decimation % 1) &&
          (problem3_DW.ToFile2_IWORK.Count*2)+1 < 100000000 ) {
        FILE *fp = (FILE *) problem3_DW.ToFile2_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[2];
          problem3_DW.ToFile2_IWORK.Decimation = 0;
          u[0] = problem3_M->Timing.t[1];
          u[1] = problem3_B.Sum4;
          if (fwrite(u, sizeof(real_T), 2, fp) != 2) {
            rtmSetErrorStatus(problem3_M,
                              "Error writing to MAT-file pitch_ref.mat");
            return;
          }

          if (((++problem3_DW.ToFile2_IWORK.Count)*2)+1 >= 100000000) {
            (void)fprintf(stdout,
                          "*** The ToFile block will stop logging data before\n"
                          "    the simulation has ended, because it has reached\n"
                          "    the maximum number of elements (100000000)\n"
                          "    allowed in MAT-file pitch_ref.mat.\n");
          }
        }
      }
    }

    /* Gain: '<S7>/Elevation: Count to rad' */
    problem3_B.ElevationCounttorad = problem3_P.ElevationCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o3;

    /* Gain: '<S10>/Gain' */
    problem3_B.Gain_e = problem3_P.Gain_Gain_lv * problem3_B.ElevationCounttorad;

    /* Sum: '<Root>/Sum' incorporates:
     *  Constant: '<Root>/elavation_offset [deg]'
     */
    problem3_B.Sum = problem3_B.Gain_e + problem3_P.elavation_offsetdeg_Value;
  }

  /* TransferFcn: '<S7>/Elevation: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += problem3_P.ElevationTransferFcn_C *
    problem3_X.ElevationTransferFcn_CSTATE;
  rtb_Backgain += problem3_P.ElevationTransferFcn_D *
    problem3_B.ElevationCounttorad;

  /* Gain: '<S11>/Gain' */
  problem3_B.Gain_dg = problem3_P.Gain_Gain_n * rtb_Backgain;

  /* Gain: '<S2>/Gain1' */
  rtb_Gain1_idx_4 = problem3_P.Gain1_Gain_f * problem3_B.Sum;
  rtb_Gain1_idx_5 = problem3_P.Gain1_Gain_f * problem3_B.Gain_dg;

  /* Sum: '<Root>/Sum1' incorporates:
   *  Constant: '<Root>/Vd_bias'
   *  Gain: '<S2>/Gain1'
   *  Gain: '<S9>/K_pd'
   *  Gain: '<S9>/K_pp'
   *  Sum: '<S9>/Sum2'
   *  Sum: '<S9>/Sum3'
   */
  problem3_B.Sum1 = ((problem3_B.Sum4 - problem3_P.Gain1_Gain_f *
                      problem3_B.Gain_i) * problem3_P.K_pp -
                     problem3_P.Gain1_Gain_f * problem3_B.Gain_b *
                     problem3_P.K_pd) + problem3_P.Vd_ff;
  if (rtmIsMajorTimeStep(problem3_M)) {
  }

  /* Integrator: '<S6>/Integrator'
   *
   * Regarding '<S6>/Integrator':
   *  Limited Integrator
   */
  if (problem3_X.Integrator_CSTATE >= problem3_P.Integrator_UpperSat ) {
    problem3_X.Integrator_CSTATE = problem3_P.Integrator_UpperSat;
  } else if (problem3_X.Integrator_CSTATE <= (problem3_P.Integrator_LowerSat) )
  {
    problem3_X.Integrator_CSTATE = (problem3_P.Integrator_LowerSat);
  }

  rtb_Backgain = problem3_X.Integrator_CSTATE;

  /* Sum: '<S6>/Sum' incorporates:
   *  Constant: '<Root>/elevation_ref'
   */
  rtb_Gain1_idx_4 = problem3_P.elevation_ref_Value - rtb_Gain1_idx_4;

  /* Sum: '<Root>/Sum2' incorporates:
   *  Constant: '<Root>/Vs_bias'
   *  Gain: '<S6>/K_ed'
   *  Gain: '<S6>/K_ep'
   *  Sum: '<S6>/Sum1'
   */
  problem3_B.Sum2 = ((problem3_P.K_ep * rtb_Gain1_idx_4 + rtb_Backgain) -
                     problem3_P.K_ed * rtb_Gain1_idx_5) + problem3_P.Vs_ff;
  if (rtmIsMajorTimeStep(problem3_M)) {
  }

  /* Gain: '<S1>/Back gain' incorporates:
   *  Sum: '<S1>/Subtract'
   */
  rtb_Backgain = (problem3_B.Sum2 - problem3_B.Sum1) * problem3_P.Backgain_Gain;

  /* Gain: '<S6>/K_ei' */
  problem3_B.K_ei = problem3_P.K_ei * rtb_Gain1_idx_4;
  if (rtmIsMajorTimeStep(problem3_M)) {
  }

  /* Derivative: '<S7>/Derivative' */
  if ((problem3_DW.TimeStampA >= problem3_M->Timing.t[0]) &&
      (problem3_DW.TimeStampB >= problem3_M->Timing.t[0])) {
    rtb_Gain1_idx_4 = 0.0;
  } else {
    rtb_Gain1_idx_4 = problem3_DW.TimeStampA;
    lastU = &problem3_DW.LastUAtTimeA;
    if (problem3_DW.TimeStampA < problem3_DW.TimeStampB) {
      if (problem3_DW.TimeStampB < problem3_M->Timing.t[0]) {
        rtb_Gain1_idx_4 = problem3_DW.TimeStampB;
        lastU = &problem3_DW.LastUAtTimeB;
      }
    } else {
      if (problem3_DW.TimeStampA >= problem3_M->Timing.t[0]) {
        rtb_Gain1_idx_4 = problem3_DW.TimeStampB;
        lastU = &problem3_DW.LastUAtTimeB;
      }
    }

    rtb_Gain1_idx_4 = (problem3_B.PitchCounttorad - *lastU) /
      (problem3_M->Timing.t[0] - rtb_Gain1_idx_4);
  }

  /* End of Derivative: '<S7>/Derivative' */

  /* Gain: '<S14>/Gain' */
  problem3_B.Gain_l = problem3_P.Gain_Gain_a1 * rtb_Gain1_idx_4;
  if (rtmIsMajorTimeStep(problem3_M)) {
  }

  /* Saturate: '<S7>/Back motor: Saturation' */
  if (rtb_Backgain > problem3_P.BackmotorSaturation_UpperSat) {
    problem3_B.BackmotorSaturation = problem3_P.BackmotorSaturation_UpperSat;
  } else if (rtb_Backgain < problem3_P.BackmotorSaturation_LowerSat) {
    problem3_B.BackmotorSaturation = problem3_P.BackmotorSaturation_LowerSat;
  } else {
    problem3_B.BackmotorSaturation = rtb_Backgain;
  }

  /* End of Saturate: '<S7>/Back motor: Saturation' */
  if (rtmIsMajorTimeStep(problem3_M)) {
  }

  /* Gain: '<S1>/Front gain' incorporates:
   *  Sum: '<S1>/Add'
   */
  rtb_Gain1_idx_4 = (problem3_B.Sum1 + problem3_B.Sum2) *
    problem3_P.Frontgain_Gain;

  /* Saturate: '<S7>/Front motor: Saturation' */
  if (rtb_Gain1_idx_4 > problem3_P.FrontmotorSaturation_UpperSat) {
    problem3_B.FrontmotorSaturation = problem3_P.FrontmotorSaturation_UpperSat;
  } else if (rtb_Gain1_idx_4 < problem3_P.FrontmotorSaturation_LowerSat) {
    problem3_B.FrontmotorSaturation = problem3_P.FrontmotorSaturation_LowerSat;
  } else {
    problem3_B.FrontmotorSaturation = rtb_Gain1_idx_4;
  }

  /* End of Saturate: '<S7>/Front motor: Saturation' */
  if (rtmIsMajorTimeStep(problem3_M)) {
    /* S-Function (hil_write_analog_block): '<S7>/HIL Write Analog' */

    /* S-Function Block: problem3/Helicopter_interface/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      problem3_DW.HILWriteAnalog_Buffer[0] = problem3_B.FrontmotorSaturation;
      problem3_DW.HILWriteAnalog_Buffer[1] = problem3_B.BackmotorSaturation;
      result = hil_write_analog(problem3_DW.HILInitialize_Card,
        problem3_P.HILWriteAnalog_channels, 2,
        &problem3_DW.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(problem3_M, _rt_error_message);
      }
    }
  }

  /* FromWorkspace: '<Root>/From Workspace' */
  {
    real_T *pDataValues = (real_T *) problem3_DW.FromWorkspace_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *) problem3_DW.FromWorkspace_PWORK.TimePtr;
    int_T currTimeIndex = problem3_DW.FromWorkspace_IWORK.PrevIndex;
    real_T t = problem3_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[140]) {
      currTimeIndex = 139;
    } else {
      if (t < pTimeValues[currTimeIndex]) {
        while (t < pTimeValues[currTimeIndex]) {
          currTimeIndex--;
        }
      } else {
        while (t >= pTimeValues[currTimeIndex + 1]) {
          currTimeIndex++;
        }
      }
    }

    problem3_DW.FromWorkspace_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          problem3_B.FromWorkspace = pDataValues[currTimeIndex];
        } else {
          problem3_B.FromWorkspace = pDataValues[currTimeIndex + 1];
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        problem3_B.FromWorkspace = (real_T) rtInterpolate(d1, d2, f1, f2);
        pDataValues += 141;
      }
    }
  }

  if (rtmIsMajorTimeStep(problem3_M)) {
    /* ToFile: '<Root>/To File3' */
    {
      if (!(++problem3_DW.ToFile3_IWORK.Decimation % 1) &&
          (problem3_DW.ToFile3_IWORK.Count*2)+1 < 100000000 ) {
        FILE *fp = (FILE *) problem3_DW.ToFile3_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[2];
          problem3_DW.ToFile3_IWORK.Decimation = 0;
          u[0] = problem3_M->Timing.t[1];
          u[1] = problem3_B.FromWorkspace;
          if (fwrite(u, sizeof(real_T), 2, fp) != 2) {
            rtmSetErrorStatus(problem3_M,
                              "Error writing to MAT-file travel_ref.mat");
            return;
          }

          if (((++problem3_DW.ToFile3_IWORK.Count)*2)+1 >= 100000000) {
            (void)fprintf(stdout,
                          "*** The ToFile block will stop logging data before\n"
                          "    the simulation has ended, because it has reached\n"
                          "    the maximum number of elements (100000000)\n"
                          "    allowed in MAT-file travel_ref.mat.\n");
          }
        }
      }
    }
  }
}

/* Model update function */
void problem3_update(void)
{
  real_T *lastU;

  /* Update for Derivative: '<S7>/Derivative' */
  if (problem3_DW.TimeStampA == (rtInf)) {
    problem3_DW.TimeStampA = problem3_M->Timing.t[0];
    lastU = &problem3_DW.LastUAtTimeA;
  } else if (problem3_DW.TimeStampB == (rtInf)) {
    problem3_DW.TimeStampB = problem3_M->Timing.t[0];
    lastU = &problem3_DW.LastUAtTimeB;
  } else if (problem3_DW.TimeStampA < problem3_DW.TimeStampB) {
    problem3_DW.TimeStampA = problem3_M->Timing.t[0];
    lastU = &problem3_DW.LastUAtTimeA;
  } else {
    problem3_DW.TimeStampB = problem3_M->Timing.t[0];
    lastU = &problem3_DW.LastUAtTimeB;
  }

  *lastU = problem3_B.PitchCounttorad;

  /* End of Update for Derivative: '<S7>/Derivative' */
  if (rtmIsMajorTimeStep(problem3_M)) {
    rt_ertODEUpdateContinuousStates(&problem3_M->solverInfo);
  }

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++problem3_M->Timing.clockTick0)) {
    ++problem3_M->Timing.clockTickH0;
  }

  problem3_M->Timing.t[0] = rtsiGetSolverStopTime(&problem3_M->solverInfo);

  {
    /* Update absolute timer for sample time: [0.002s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick1"
     * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++problem3_M->Timing.clockTick1)) {
      ++problem3_M->Timing.clockTickH1;
    }

    problem3_M->Timing.t[1] = problem3_M->Timing.clockTick1 *
      problem3_M->Timing.stepSize1 + problem3_M->Timing.clockTickH1 *
      problem3_M->Timing.stepSize1 * 4294967296.0;
  }
}

/* Derivatives for root system: '<Root>' */
void problem3_derivatives(void)
{
  XDot_problem3_T *_rtXdot;
  _rtXdot = ((XDot_problem3_T *) problem3_M->ModelData.derivs);

  /* Derivatives for TransferFcn: '<S7>/Travel: Transfer Fcn' */
  _rtXdot->TravelTransferFcn_CSTATE = 0.0;
  _rtXdot->TravelTransferFcn_CSTATE += problem3_P.TravelTransferFcn_A *
    problem3_X.TravelTransferFcn_CSTATE;
  _rtXdot->TravelTransferFcn_CSTATE += problem3_B.TravelCounttorad;

  /* Derivatives for TransferFcn: '<S7>/Pitch: Transfer Fcn' */
  _rtXdot->PitchTransferFcn_CSTATE = 0.0;
  _rtXdot->PitchTransferFcn_CSTATE += problem3_P.PitchTransferFcn_A *
    problem3_X.PitchTransferFcn_CSTATE;
  _rtXdot->PitchTransferFcn_CSTATE += problem3_B.PitchCounttorad;

  /* Derivatives for TransferFcn: '<S7>/Elevation: Transfer Fcn' */
  _rtXdot->ElevationTransferFcn_CSTATE = 0.0;
  _rtXdot->ElevationTransferFcn_CSTATE += problem3_P.ElevationTransferFcn_A *
    problem3_X.ElevationTransferFcn_CSTATE;
  _rtXdot->ElevationTransferFcn_CSTATE += problem3_B.ElevationCounttorad;

  /* Derivatives for Integrator: '<S6>/Integrator' */
  {
    boolean_T lsat;
    boolean_T usat;
    lsat = ( problem3_X.Integrator_CSTATE <= (problem3_P.Integrator_LowerSat) );
    usat = ( problem3_X.Integrator_CSTATE >= problem3_P.Integrator_UpperSat );
    if ((!lsat && !usat) ||
        (lsat && (problem3_B.K_ei > 0)) ||
        (usat && (problem3_B.K_ei < 0)) ) {
      ((XDot_problem3_T *) problem3_M->ModelData.derivs)->Integrator_CSTATE =
        problem3_B.K_ei;
    } else {
      /* in saturation */
      ((XDot_problem3_T *) problem3_M->ModelData.derivs)->Integrator_CSTATE =
        0.0;
    }
  }
}

/* Model initialize function */
void problem3_initialize(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: problem3/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q8_usb", "0", &problem3_DW.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(problem3_M, _rt_error_message);
      return;
    }

    is_switching = false;
    result = hil_set_card_specific_options(problem3_DW.HILInitialize_Card,
      "update_rate=normal;decimation=1", 32);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(problem3_M, _rt_error_message);
      return;
    }

    result = hil_watchdog_clear(problem3_DW.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(problem3_M, _rt_error_message);
      return;
    }

    if ((problem3_P.HILInitialize_set_analog_input_ && !is_switching) ||
        (problem3_P.HILInitialize_set_analog_inpu_m && is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums = &problem3_DW.HILInitialize_AIMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMinimums[i1] = problem3_P.HILInitialize_analog_input_mini;
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums = &problem3_DW.HILInitialize_AIMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMaximums[i1] = problem3_P.HILInitialize_analog_input_maxi;
        }
      }

      result = hil_set_analog_input_ranges(problem3_DW.HILInitialize_Card,
        problem3_P.HILInitialize_analog_input_chan, 8U,
        &problem3_DW.HILInitialize_AIMinimums[0],
        &problem3_DW.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(problem3_M, _rt_error_message);
        return;
      }
    }

    if ((problem3_P.HILInitialize_set_analog_output && !is_switching) ||
        (problem3_P.HILInitialize_set_analog_outp_b && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOMinimums = &problem3_DW.HILInitialize_AOMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMinimums[i1] = problem3_P.HILInitialize_analog_output_min;
        }
      }

      {
        int_T i1;
        real_T *dw_AOMaximums = &problem3_DW.HILInitialize_AOMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMaximums[i1] = problem3_P.HILInitialize_analog_output_max;
        }
      }

      result = hil_set_analog_output_ranges(problem3_DW.HILInitialize_Card,
        problem3_P.HILInitialize_analog_output_cha, 8U,
        &problem3_DW.HILInitialize_AOMinimums[0],
        &problem3_DW.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(problem3_M, _rt_error_message);
        return;
      }
    }

    if ((problem3_P.HILInitialize_set_analog_outp_e && !is_switching) ||
        (problem3_P.HILInitialize_set_analog_outp_j && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &problem3_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = problem3_P.HILInitialize_initial_analog_ou;
        }
      }

      result = hil_write_analog(problem3_DW.HILInitialize_Card,
        problem3_P.HILInitialize_analog_output_cha, 8U,
        &problem3_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(problem3_M, _rt_error_message);
        return;
      }
    }

    if (problem3_P.HILInitialize_set_analog_outp_p) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &problem3_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = problem3_P.HILInitialize_watchdog_analog_o;
        }
      }

      result = hil_watchdog_set_analog_expiration_state
        (problem3_DW.HILInitialize_Card,
         problem3_P.HILInitialize_analog_output_cha, 8U,
         &problem3_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(problem3_M, _rt_error_message);
        return;
      }
    }

    if ((problem3_P.HILInitialize_set_encoder_param && !is_switching) ||
        (problem3_P.HILInitialize_set_encoder_par_m && is_switching)) {
      {
        int_T i1;
        int32_T *dw_QuadratureModes =
          &problem3_DW.HILInitialize_QuadratureModes[0];
        for (i1=0; i1 < 8; i1++) {
          dw_QuadratureModes[i1] = problem3_P.HILInitialize_quadrature;
        }
      }

      result = hil_set_encoder_quadrature_mode(problem3_DW.HILInitialize_Card,
        problem3_P.HILInitialize_encoder_channels, 8U,
        (t_encoder_quadrature_mode *)
        &problem3_DW.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(problem3_M, _rt_error_message);
        return;
      }
    }

    if ((problem3_P.HILInitialize_set_encoder_count && !is_switching) ||
        (problem3_P.HILInitialize_set_encoder_cou_k && is_switching)) {
      {
        int_T i1;
        int32_T *dw_InitialEICounts =
          &problem3_DW.HILInitialize_InitialEICounts[0];
        for (i1=0; i1 < 8; i1++) {
          dw_InitialEICounts[i1] = problem3_P.HILInitialize_initial_encoder_c;
        }
      }

      result = hil_set_encoder_counts(problem3_DW.HILInitialize_Card,
        problem3_P.HILInitialize_encoder_channels, 8U,
        &problem3_DW.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(problem3_M, _rt_error_message);
        return;
      }
    }

    if ((problem3_P.HILInitialize_set_pwm_params_at && !is_switching) ||
        (problem3_P.HILInitialize_set_pwm_params__f && is_switching)) {
      uint32_T num_duty_cycle_modes = 0;
      uint32_T num_frequency_modes = 0;

      {
        int_T i1;
        int32_T *dw_POModeValues = &problem3_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = problem3_P.HILInitialize_pwm_modes;
        }
      }

      result = hil_set_pwm_mode(problem3_DW.HILInitialize_Card,
        problem3_P.HILInitialize_pwm_channels, 8U, (t_pwm_mode *)
        &problem3_DW.HILInitialize_POModeValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(problem3_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        const uint32_T *p_HILInitialize_pwm_channels =
          problem3_P.HILInitialize_pwm_channels;
        int32_T *dw_POModeValues = &problem3_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          if (dw_POModeValues[i1] == PWM_DUTY_CYCLE_MODE || dw_POModeValues[i1] ==
              PWM_ONE_SHOT_MODE || dw_POModeValues[i1] == PWM_TIME_MODE) {
            problem3_DW.HILInitialize_POSortedChans[num_duty_cycle_modes] =
              p_HILInitialize_pwm_channels[i1];
            problem3_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes] =
              problem3_P.HILInitialize_pwm_frequency;
            num_duty_cycle_modes++;
          } else {
            problem3_DW.HILInitialize_POSortedChans[7U - num_frequency_modes] =
              p_HILInitialize_pwm_channels[i1];
            problem3_DW.HILInitialize_POSortedFreqs[7U - num_frequency_modes] =
              problem3_P.HILInitialize_pwm_frequency;
            num_frequency_modes++;
          }
        }
      }

      if (num_duty_cycle_modes > 0) {
        result = hil_set_pwm_frequency(problem3_DW.HILInitialize_Card,
          &problem3_DW.HILInitialize_POSortedChans[0], num_duty_cycle_modes,
          &problem3_DW.HILInitialize_POSortedFreqs[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(problem3_M, _rt_error_message);
          return;
        }
      }

      if (num_frequency_modes > 0) {
        result = hil_set_pwm_duty_cycle(problem3_DW.HILInitialize_Card,
          &problem3_DW.HILInitialize_POSortedChans[num_duty_cycle_modes],
          num_frequency_modes,
          &problem3_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(problem3_M, _rt_error_message);
          return;
        }
      }

      {
        int_T i1;
        int32_T *dw_POModeValues = &problem3_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = problem3_P.HILInitialize_pwm_configuration;
        }
      }

      {
        int_T i1;
        int32_T *dw_POAlignValues = &problem3_DW.HILInitialize_POAlignValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POAlignValues[i1] = problem3_P.HILInitialize_pwm_alignment;
        }
      }

      {
        int_T i1;
        int32_T *dw_POPolarityVals = &problem3_DW.HILInitialize_POPolarityVals[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POPolarityVals[i1] = problem3_P.HILInitialize_pwm_polarity;
        }
      }

      result = hil_set_pwm_configuration(problem3_DW.HILInitialize_Card,
        problem3_P.HILInitialize_pwm_channels, 8U,
        (t_pwm_configuration *) &problem3_DW.HILInitialize_POModeValues[0],
        (t_pwm_alignment *) &problem3_DW.HILInitialize_POAlignValues[0],
        (t_pwm_polarity *) &problem3_DW.HILInitialize_POPolarityVals[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(problem3_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        real_T *dw_POSortedFreqs = &problem3_DW.HILInitialize_POSortedFreqs[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POSortedFreqs[i1] = problem3_P.HILInitialize_pwm_leading_deadb;
        }
      }

      {
        int_T i1;
        real_T *dw_POValues = &problem3_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = problem3_P.HILInitialize_pwm_trailing_dead;
        }
      }

      result = hil_set_pwm_deadband(problem3_DW.HILInitialize_Card,
        problem3_P.HILInitialize_pwm_channels, 8U,
        &problem3_DW.HILInitialize_POSortedFreqs[0],
        &problem3_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(problem3_M, _rt_error_message);
        return;
      }
    }

    if ((problem3_P.HILInitialize_set_pwm_outputs_a && !is_switching) ||
        (problem3_P.HILInitialize_set_pwm_outputs_g && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &problem3_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = problem3_P.HILInitialize_initial_pwm_outpu;
        }
      }

      result = hil_write_pwm(problem3_DW.HILInitialize_Card,
        problem3_P.HILInitialize_pwm_channels, 8U,
        &problem3_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(problem3_M, _rt_error_message);
        return;
      }
    }

    if (problem3_P.HILInitialize_set_pwm_outputs_o) {
      {
        int_T i1;
        real_T *dw_POValues = &problem3_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = problem3_P.HILInitialize_watchdog_pwm_outp;
        }
      }

      result = hil_watchdog_set_pwm_expiration_state
        (problem3_DW.HILInitialize_Card, problem3_P.HILInitialize_pwm_channels,
         8U, &problem3_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(problem3_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for S-Function (hil_read_encoder_timebase_block): '<S7>/HIL Read Encoder Timebase' */

  /* S-Function Block: problem3/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_create_encoder_reader(problem3_DW.HILInitialize_Card,
      problem3_P.HILReadEncoderTimebase_samples_,
      problem3_P.HILReadEncoderTimebase_channels, 3,
      &problem3_DW.HILReadEncoderTimebase_Task);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(problem3_M, _rt_error_message);
    }
  }

  /* Start for ToFile: '<Root>/To File' */
  {
    char fileName[509] = "pitch.mat";
    FILE *fp = (NULL);
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(problem3_M, "Error creating .mat file pitch.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp,2,0,"measured_pitch")) {
      rtmSetErrorStatus(problem3_M,
                        "Error writing mat file header to file pitch.mat");
      return;
    }

    problem3_DW.ToFile_IWORK.Count = 0;
    problem3_DW.ToFile_IWORK.Decimation = -1;
    problem3_DW.ToFile_PWORK.FilePtr = fp;
  }

  /* Start for ToFile: '<Root>/To File1' */
  {
    char fileName[509] = "travel.mat";
    FILE *fp = (NULL);
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(problem3_M, "Error creating .mat file travel.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp,2,0,"measured_travel")) {
      rtmSetErrorStatus(problem3_M,
                        "Error writing mat file header to file travel.mat");
      return;
    }

    problem3_DW.ToFile1_IWORK.Count = 0;
    problem3_DW.ToFile1_IWORK.Decimation = -1;
    problem3_DW.ToFile1_PWORK.FilePtr = fp;
  }

  /* Start for FromWorkspace: '<S8>/From Workspace3' */
  {
    static real_T pTimeValues0[] = { 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75,
      3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0,
      6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25, 8.5, 8.75, 9.0, 9.25,
      9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25, 11.5, 11.75, 12.0, 12.25,
      12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0, 14.25, 14.5, 14.75, 15.0,
      15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75, 17.0, 17.25, 17.5, 17.75,
      18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5, 19.75, 20.0, 20.25, 20.5,
      20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25, 22.5, 22.75, 23.0, 23.25,
      23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0, 25.25, 25.5, 25.75, 26.0,
      26.25, 26.5, 26.75, 27.0, 27.25, 27.5, 27.75, 28.0, 28.25, 28.5, 28.75,
      29.0, 29.25, 29.5, 29.75, 30.0, 30.25, 30.5, 30.75, 31.0, 31.25, 31.5,
      31.75, 32.0, 32.25, 32.5, 32.75, 33.0, 33.25, 33.5, 33.75, 34.0, 34.25,
      34.5, 34.75, 35.0, 35.25, 35.5, 35.75, 36.0 } ;

    static real_T pDataValues0[] = { 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.137842141357523, 3.1262155534531373,
      3.1033093000254177, 3.0666274151872006, 3.0144539223910845,
      2.9456562771158223, 2.8595077632936725, 2.7555515879678714,
      2.6335051104965794, 2.4931956060433924, 2.3345185760823788,
      2.1574113214985227, 1.96183648407086, 1.7516021122598986,
      1.5345792086565013, 1.3182739414574947, 1.108562117396344,
      0.90963650800063778, 0.72433173041414511, 0.55449755010832269,
      0.40131123741616709, 0.26550823564839071, 0.14754140816898981,
      0.047686433076618941, -0.033890434885039504, -0.097089708132804431,
      -0.14185224237524749, -0.1681431789178662, -0.17739386917910949,
      -0.172487569853685, -0.15710545639363724, -0.13505545820928427,
      -0.10978472456783928, -0.08410200575622831, -0.060073243705058232,
      -0.039038680112135675, -0.021703755768764897, -0.0082661385761011938,
      0.0014480180521598508, 0.0078560697682221771, 0.011509140048521877,
      0.013005114270146353, 0.012924992403316005, 0.011790760213555082,
      0.010041552293830643, 0.0080242911898826784, 0.0059949552401125155,
      0.0041269557824722674, 0.002523635906620587, 0.001232520268376065,
      0.00025956499368645424, -0.00041777645181630639, -0.00083937289969922637,
      -0.0010534372648807737, -0.00110978078857205, -0.0010551342373548865,
      -0.00093025223982172315, -0.00076846079828918367, -0.00059529849725052631,
      -0.00042892465495781565, -0.00028101023259235943, -0.00015787971043845165,
      -6.1726788911796064E-5, 8.2216193074516516E-6, 5.4673385948917387E-5,
      8.1414237098516944E-5, 9.2640319787744009E-5, 9.2474170510006867E-5,
      8.46446726442204E-5, 7.2304234787592787E-5, 5.7953699453178867E-5,
      4.3446117997920216E-5, 3.0043358600500043E-5, 1.8503594410534691E-5,
      9.1823006751724239E-6, 2.1339268228718973E-6, -2.794462010272638E-6,
      -5.8827709807183229E-6, -7.47308808621365E-6, -7.9208366974172172E-6,
      -7.5607499912566815E-6, -6.685677866308858E-6, -5.53581569708402E-6,
      -4.2958614006927632E-6, -3.09775914344273E-6, -2.0269851623538277E-6,
      -1.1307011921095014E-6, -4.2648898386397463E-7, 8.9253070272431473E-8,
      4.3384723708263285E-7, 6.3213765142667077E-7, 7.11679557472354E-7,
      6.9929177222740138E-7, 6.1884802058298051E-7, 4.9012269740324077E-7,
      3.2847518897083036E-7, 1.4514398304184174E-7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.01500204890654452,
      -0.046506351615187548, -0.091625013708521649, -0.14672753935051416,
      -0.20869397118210767, -0.27519058109869365, -0.34459405528624421,
      -0.41582470130085009, -0.48818590988281235, -0.56123801781039206,
      -0.63470811984169739, -0.70842901833307059, -0.78229934970829484,
      -0.84093748724148887, -0.8680916144112335, -0.86522106879367122,
      -0.83884729624224619, -0.79570243758046988, -0.74121911034361487,
      -0.679336721220934, -0.61274525076626674, -0.54321200706874984,
      -0.471867309915248, -0.39941990036712782, -0.32630747184427811,
      -0.25279709298870406, -0.17905013696741659, -0.10516374616811908,
      -0.037002761042617474, 0.019625197304053647, 0.061528453842546696,
      0.088199992739767516, 0.10108293456813565, 0.10273087524879962,
      0.09611504820703598, 0.084138254374045915, 0.0693396973758388,
      0.053750468773010505, 0.038856626515399867, 0.025632206866604997,
      0.014612281123554489, 0.0059838968888535982, -0.00032048746496569691,
      -0.0045369287566880074, -0.0069968316765420631, -0.0080690444134361612,
      -0.00811734379672496, -0.0074719978282053014, -0.0064132795010510332,
      -0.0051644625506223956, -0.0038918210964027523, -0.0027093657796553516,
      -0.0016863857891759885, -0.00085625745837049848, -0.00022537409240941348,
      0.00021858620722434532, 0.0004995279924883447, 0.00064716576848584929,
      0.00069264920651032061, 0.00066549537152653393, 0.00059165769181751613,
      0.00049252209097132228, 0.00038461168846231361, 0.00027979363523268207,
      0.00018580706892155422, 0.00010696340695408943, 4.4904333112599546E-5,
      -6.6459475525732754E-7, -3.1317989107454563E-5, -4.9361749070819279E-5,
      -5.7402138981964409E-5, -5.8030323465343369E-5, -5.3611035233989458E-5,
      -4.6159054404170164E-5, -3.7285172585757818E-5, -2.8193493053510857E-5,
      -1.9713552976886895E-5, -1.2353233526091492E-5, -6.3612660662900613E-6,
      -1.7909920891230261E-6, 1.4403491803333944E-6, 3.5002908554825408E-6,
      4.5994510325906022E-6, 4.95981954125627E-6, 4.7924113846913812E-6,
      4.2830982800468573E-6, 3.5851382366685512E-6, 2.8168511886733544E-6,
      2.0629705722368719E-6, 1.3783790229320528E-6, 7.9316401306739872E-7,
      3.1816997987398048E-7, -4.9548785288563654E-8, -3.2177265088643595E-7,
      -5.1489893702771152E-7, -6.4658767803839438E-7, -7.3332246802470733E-7,
      -7.8889491054332731E-7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.10602875205758999, 0.22266037932046823, 0.31888147181165888,
      0.38944360630438574, 0.43795507376741688, 0.46997264229199831,
      0.4905172487602496, 0.50343100139455166, 0.51142138583254448,
      0.51630439853687438, 0.51925862122035837, 0.5210311548151294,
      0.52208728936630211, 0.41443196081916378, 0.1919150000445301,
      -0.020287920096817472, -0.18639975163753869, -0.30493138299849926,
      -0.38506734846286844, -0.4373610919984362, -0.47064308034546981,
      -0.49143441009445715, -0.50423707128524442, -0.5120306213389525,
      -0.51673072146616084, -0.51954328243138836, -0.52121531948434141,
      -0.52220079124193863, -0.48173581060484882, -0.40022478206265588,
      -0.296156213387799, -0.18850425044574579, -0.091051712548316066,
      -0.011647015344018965, 0.046758138793423977, 0.08464740459617949,
      0.1045905489510945, 0.11017871387672144, 0.10526398877290666,
      0.0934651473650344, 0.0778846263866169, 0.060982124390805824,
      0.044556980822223695, 0.029800196693233791, 0.017385654343638243,
      0.0075779901215631135, 0.00034136159455317925, -0.0045610588327712771,
      -0.0074826167871035569, -0.0088261612533308432, -0.0089945437469947549,
      -0.0083571425715198572, -0.0072300318727855058, -0.0058670299967350289,
      -0.0044588426814908483, -0.0031377418389463749, -0.001985589239934824,
      -0.0010434474143340446, -0.00032145956865307635, 0.00019191293510899219,
      0.000521856519983204, 0.00070065256471757789, 0.00076266951156995632,
      0.000740813967903274, 0.00066426115515951737, 0.00055723689068510631,
      0.00043860983223445124, 0.00032206377846225737, 0.00021664648416197857,
      0.00012752640421372292, 5.6826405203812169E-5, 4.4397681336333263E-6,
      -3.1233841267129036E-5, -5.2667754297397535E-5, -6.2717207427786674E-5,
      -6.4256518487168268E-5, -5.9932977652320057E-5, -5.2019926678064509E-5,
      -4.2348937411100649E-5, -3.2300950853575868E-5, -2.2837885844951756E-5,
      -1.4558880897078396E-5, -7.76844421367555E-6, -2.5469468628084791E-6,
      1.1831770368841962E-6, 3.5996306284766E-6, 4.9329151428030227E-6,
      5.4299595519200184E-6, 5.32814039581358E-6, 4.8384317324297511E-6,
      4.1360763058733689E-6, 3.357076538710032E-6, 2.5988959151822776E-6,
      1.923974436777529E-6, 1.3649429543795644E-6, 9.30725849820839E-7,
      6.1300849293776181E-7, 3.9276479931610446E-7, 2.4660623571375195E-7, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.42411500823222226,
      0.46652650905386861, 0.3848843699671185, 0.282248537973263,
      0.19404586985448027, 0.12807027410068142, 0.082178425875360783,
      0.051655010539563805, 0.031961537754327259, 0.019532050819675077,
      0.0118168907362918, 0.0070901343814398941, 0.0042245382070466133,
      -0.43062131418619759, -0.89006784309617892, -0.84881168056303469,
      -0.66444732616052926, -0.47412652544148659, -0.32054386185512113,
      -0.2091749741399154, -0.13312795338577868, -0.0831653189935938,
      -0.051210644760793143, -0.031174200212476862, -0.018800400506477324,
      -0.011250243858554609, -0.00668814820945643, -0.0039418870280330471,
      0.16185992255071494, 0.32604411417112728, 0.41627427470178319,
      0.43060785177056865, 0.38981015159207455, 0.31761878881954408,
      0.23362061655212743, 0.15155706321337781, 0.079772577422015739,
      0.02235265970486346, -0.019658900412903427, -0.047195365629133332,
      -0.062322083911314309, -0.06761000798088862, -0.065700574271972831,
      -0.059027136513603924, -0.049658169396026505, -0.039230656885944827,
      -0.028946514105684048, -0.019609681706942134, -0.011686231814973426,
      -0.005374177862553454, -0.00067352997229995363, 0.0025496047042552874,
      0.0045084427972930958, 0.0054520075065576, 0.0056327492633324169,
      0.005284403372533583, 0.0046086103984018939, 0.0037685673047588085,
      0.0028879513850795647, 0.0020534900174039653, 0.0013197743418525382,
      0.00071518418129318717, 0.000248067789765205, -8.7422172311037838E-5,
      -0.00030621124861933532, -0.00042809705554195293, -0.00047450823144692917,
      -0.00046618421273308421, -0.00042166917484542384, -0.00035648031743733147,
      -0.00028279999368395172, -0.00020954654592502412, -0.00014269443524735821,
      -8.5735649765382739E-5, -4.01978101658653E-5, -6.1572418818351727E-6,
      1.7294165695084105E-5, 3.165220625271345E-5, 3.8683959423546688E-5,
      4.0191948585790363E-5, 3.7852262390187705E-5, 3.3116022147184686E-5,
      2.7161749089302625E-5, 2.0885991759159535E-5, 1.4920497954461949E-5,
      9.6658167220608612E-6, 5.3331404129969386E-6, 1.9881799921592343E-6,
      -4.0727426873450779E-7, -1.9588322978440704E-6, -2.8094193505342807E-6,
      -3.1159967129621E-6, -3.0327201384197715E-6, -2.6996835579277464E-6,
      -2.2361235739006107E-6, -1.7368660625436545E-6, -1.2708670718410612E-6,
      -8.80972418795382E-7, -5.8463189871816278E-7, -3.7660267731786097E-7, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0 } ;

    problem3_DW.FromWorkspace3_PWORK.TimePtr = (void *) pTimeValues0;
    problem3_DW.FromWorkspace3_PWORK.DataPtr = (void *) pDataValues0;
    problem3_DW.FromWorkspace3_IWORK.PrevIndex = 0;
  }

  /* Start for FromWorkspace: '<S8>/From Workspace4' */
  {
    static real_T pTimeValues0[] = { 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75,
      3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0,
      6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25, 8.5, 8.75, 9.0, 9.25,
      9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25, 11.5, 11.75, 12.0, 12.25,
      12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0, 14.25, 14.5, 14.75, 15.0,
      15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75, 17.0, 17.25, 17.5, 17.75,
      18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5, 19.75, 20.0, 20.25, 20.5,
      20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25, 22.5, 22.75, 23.0, 23.25,
      23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0, 25.25, 25.5, 25.75, 26.0,
      26.25, 26.5, 26.75, 27.0, 27.25, 27.5, 27.75, 28.0, 28.25, 28.5, 28.75,
      29.0, 29.25, 29.5, 29.75, 30.0, 30.25, 30.5, 30.75, 31.0, 31.25, 31.5,
      31.75, 32.0, 32.25, 32.5, 32.75, 33.0, 33.25, 33.5, 33.75, 34.0, 34.25,
      34.5, 34.75, 35.0, 35.25, 35.5, 35.75, 36.0 } ;

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.52359877559255108,
      0.5235987755908964, 0.52359877558797452, 0.52359877558519785,
      0.52359877558330559, 0.52359877557983414, 0.52359877557380263,
      0.5235987755634891, 0.52359877554821355, 0.52359877552062839,
      0.52359877552662382, 0.52359877545748013, 0.52359877525342324,
      -0.011121645196009203, -0.523598774435771, -0.52359877529878274,
      -0.523598775395403, -0.5235987754394763, -0.52359877548086287,
      -0.523598775537417, -0.523598775591149, -0.52359877552576473,
      -0.52359877523874621, -0.52359877520373821, -0.52359877559061452,
      -0.52359877419083767, -0.52359877557198931, -0.52358411231929891,
      -0.32090160189002276, -0.13965989999016623, -0.0080693127251571413,
      0.079997963984591219, 0.13192942415898806, 0.15549300098415109,
      0.15815660556039118, 0.14661841867578956, 0.12653205372991511,
      0.10239481442400881, 0.077560713907255135, 0.054339855123680565,
      0.034149732784609005, 0.017690086140494537, 0.0051197208606859361,
      -0.0037797016444124005, -0.0094621003466481956, -0.012502074848077676,
      -0.013507491852916919, -0.013057849753120651, -0.011665136734334886,
      -0.0097531016515073673, -0.0076506565255785326, -0.0055953492788059535,
      -0.0037433310003929716, -0.0021828669854603911, -0.00094910778567867406,
      -3.8476362068092564E-5, 0.0005784043049639827, 0.00094584613981720119,
      0.0011145276149831522, 0.001135188236006556, 0.0010543741819543551,
      0.000911921262254438, 0.00073981870457056709, 0.00056209830125130583,
      0.00039542379118927279, 0.00025010293975495136, 0.00013130000946463352,
      4.0282086106694737E-5, -2.4415792076953013E-5, -6.5977458664793451E-5,
      -8.8479397393997555E-5, -9.6259458350167256E-5, -9.3469622913444508E-5,
      -8.3790118731513617E-5, -7.02761218960786E-5, -5.5306597489213293E-5,
      -4.0606234880693458E-5, -2.7314805446325172E-5, -1.6082685457069173E-5,
      -7.1760371431094953E-6, -5.7972085370986409E-7, 3.9099056748732483E-6,
      6.6067373087050711E-6, 7.87298914009092E-6, 8.0734035127736711E-6,
      7.5441207995279938E-6, 6.5739876796850971E-6, 5.3957588145721775E-6,
      4.1846504661702405E-6, 3.0619262295833008E-6, 2.1015499749839566E-6,
      1.3383614879472896E-6, 7.7666570433548806E-7, 3.98541150799055E-7,
      1.715411766267733E-7, 5.5757368535904371E-8, 1.0396552269912437E-8,
      1.3163145827022966E-13, 1.26319555547681E-13, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0 } ;

    problem3_DW.FromWorkspace4_PWORK.TimePtr = (void *) pTimeValues0;
    problem3_DW.FromWorkspace4_PWORK.DataPtr = (void *) pDataValues0;
    problem3_DW.FromWorkspace4_IWORK.PrevIndex = 0;
  }

  /* Start for ToFile: '<Root>/To File2' */
  {
    char fileName[509] = "pitch_ref.mat";
    FILE *fp = (NULL);
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(problem3_M, "Error creating .mat file pitch_ref.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp,2,0,"pitch_ref")) {
      rtmSetErrorStatus(problem3_M,
                        "Error writing mat file header to file pitch_ref.mat");
      return;
    }

    problem3_DW.ToFile2_IWORK.Count = 0;
    problem3_DW.ToFile2_IWORK.Decimation = -1;
    problem3_DW.ToFile2_PWORK.FilePtr = fp;
  }

  /* Start for FromWorkspace: '<Root>/From Workspace' */
  {
    static real_T pTimeValues0[] = { 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75,
      3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0,
      6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25, 8.5, 8.75, 9.0, 9.25,
      9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25, 11.5, 11.75, 12.0, 12.25,
      12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0, 14.25, 14.5, 14.75, 15.0,
      15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75, 17.0, 17.25, 17.5, 17.75,
      18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5, 19.75, 20.0, 20.25, 20.5,
      20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25, 22.5, 22.75, 23.0, 23.25,
      23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0, 25.25, 25.5, 25.75, 26.0,
      26.25, 26.5, 26.75, 27.0, 27.25, 27.5, 27.75, 28.0, 28.25, 28.5, 28.75,
      29.0, 29.25, 29.5, 29.75, 30.0, 30.25, 30.5, 30.75, 31.0, 31.25, 31.5,
      31.75, 32.0, 32.25, 32.5, 32.75, 33.0, 33.25, 33.5, 33.75, 34.0, 34.25,
      34.5, 34.75, 35.0, 35.25, 35.5, 35.75, 36.0 } ;

    static real_T pDataValues0[] = { 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.137842141357523, 3.1262155534531373,
      3.1033093000254177, 3.0666274151872006, 3.0144539223910845,
      2.9456562771158223, 2.8595077632936725, 2.7555515879678714,
      2.6335051104965794, 2.4931956060433924, 2.3345185760823788,
      2.1574113214985227, 1.96183648407086, 1.7516021122598986,
      1.5345792086565013, 1.3182739414574947, 1.108562117396344,
      0.90963650800063778, 0.72433173041414511, 0.55449755010832269,
      0.40131123741616709, 0.26550823564839071, 0.14754140816898981,
      0.047686433076618941, -0.033890434885039504, -0.097089708132804431,
      -0.14185224237524749, -0.1681431789178662, -0.17739386917910949,
      -0.172487569853685, -0.15710545639363724, -0.13505545820928427,
      -0.10978472456783928, -0.08410200575622831, -0.060073243705058232,
      -0.039038680112135675, -0.021703755768764897, -0.0082661385761011938,
      0.0014480180521598508, 0.0078560697682221771, 0.011509140048521877,
      0.013005114270146353, 0.012924992403316005, 0.011790760213555082,
      0.010041552293830643, 0.0080242911898826784, 0.0059949552401125155,
      0.0041269557824722674, 0.002523635906620587, 0.001232520268376065,
      0.00025956499368645424, -0.00041777645181630639, -0.00083937289969922637,
      -0.0010534372648807737, -0.00110978078857205, -0.0010551342373548865,
      -0.00093025223982172315, -0.00076846079828918367, -0.00059529849725052631,
      -0.00042892465495781565, -0.00028101023259235943, -0.00015787971043845165,
      -6.1726788911796064E-5, 8.2216193074516516E-6, 5.4673385948917387E-5,
      8.1414237098516944E-5, 9.2640319787744009E-5, 9.2474170510006867E-5,
      8.46446726442204E-5, 7.2304234787592787E-5, 5.7953699453178867E-5,
      4.3446117997920216E-5, 3.0043358600500043E-5, 1.8503594410534691E-5,
      9.1823006751724239E-6, 2.1339268228718973E-6, -2.794462010272638E-6,
      -5.8827709807183229E-6, -7.47308808621365E-6, -7.9208366974172172E-6,
      -7.5607499912566815E-6, -6.685677866308858E-6, -5.53581569708402E-6,
      -4.2958614006927632E-6, -3.09775914344273E-6, -2.0269851623538277E-6,
      -1.1307011921095014E-6, -4.2648898386397463E-7, 8.9253070272431473E-8,
      4.3384723708263285E-7, 6.3213765142667077E-7, 7.11679557472354E-7,
      6.9929177222740138E-7, 6.1884802058298051E-7, 4.9012269740324077E-7,
      3.2847518897083036E-7, 1.4514398304184174E-7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

    problem3_DW.FromWorkspace_PWORK.TimePtr = (void *) pTimeValues0;
    problem3_DW.FromWorkspace_PWORK.DataPtr = (void *) pDataValues0;
    problem3_DW.FromWorkspace_IWORK.PrevIndex = 0;
  }

  /* Start for ToFile: '<Root>/To File3' */
  {
    char fileName[509] = "travel_ref.mat";
    FILE *fp = (NULL);
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(problem3_M, "Error creating .mat file travel_ref.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp,2,0,"travel_ref")) {
      rtmSetErrorStatus(problem3_M,
                        "Error writing mat file header to file travel_ref.mat");
      return;
    }

    problem3_DW.ToFile3_IWORK.Count = 0;
    problem3_DW.ToFile3_IWORK.Decimation = -1;
    problem3_DW.ToFile3_PWORK.FilePtr = fp;
  }

  /* InitializeConditions for TransferFcn: '<S7>/Travel: Transfer Fcn' */
  problem3_X.TravelTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S7>/Pitch: Transfer Fcn' */
  problem3_X.PitchTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S7>/Elevation: Transfer Fcn' */
  problem3_X.ElevationTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S6>/Integrator' */
  problem3_X.Integrator_CSTATE = problem3_P.Integrator_IC;

  /* InitializeConditions for Derivative: '<S7>/Derivative' */
  problem3_DW.TimeStampA = (rtInf);
  problem3_DW.TimeStampB = (rtInf);
}

/* Model terminate function */
void problem3_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: problem3/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_pwm_outputs = 0;
    hil_task_stop_all(problem3_DW.HILInitialize_Card);
    hil_monitor_stop_all(problem3_DW.HILInitialize_Card);
    is_switching = false;
    if ((problem3_P.HILInitialize_set_analog_out_ex && !is_switching) ||
        (problem3_P.HILInitialize_set_analog_outp_c && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &problem3_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = problem3_P.HILInitialize_final_analog_outp;
        }
      }

      num_final_analog_outputs = 8U;
    }

    if ((problem3_P.HILInitialize_set_pwm_output_ap && !is_switching) ||
        (problem3_P.HILInitialize_set_pwm_outputs_p && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &problem3_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = problem3_P.HILInitialize_final_pwm_outputs;
        }
      }

      num_final_pwm_outputs = 8U;
    }

    if (0
        || num_final_analog_outputs > 0
        || num_final_pwm_outputs > 0
        ) {
      /* Attempt to write the final outputs atomically (due to firmware issue in old Q2-USB). Otherwise write channels individually */
      result = hil_write(problem3_DW.HILInitialize_Card
                         , problem3_P.HILInitialize_analog_output_cha,
                         num_final_analog_outputs
                         , problem3_P.HILInitialize_pwm_channels,
                         num_final_pwm_outputs
                         , NULL, 0
                         , NULL, 0
                         , &problem3_DW.HILInitialize_AOVoltages[0]
                         , &problem3_DW.HILInitialize_POValues[0]
                         , (t_boolean *) NULL
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog(problem3_DW.HILInitialize_Card,
            problem3_P.HILInitialize_analog_output_cha, num_final_analog_outputs,
            &problem3_DW.HILInitialize_AOVoltages[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_pwm_outputs > 0) {
          local_result = hil_write_pwm(problem3_DW.HILInitialize_Card,
            problem3_P.HILInitialize_pwm_channels, num_final_pwm_outputs,
            &problem3_DW.HILInitialize_POValues[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(problem3_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(problem3_DW.HILInitialize_Card);
    hil_monitor_delete_all(problem3_DW.HILInitialize_Card);
    hil_close(problem3_DW.HILInitialize_Card);
    problem3_DW.HILInitialize_Card = NULL;
  }

  /* Terminate for ToFile: '<Root>/To File' */
  {
    FILE *fp = (FILE *) problem3_DW.ToFile_PWORK.FilePtr;
    if (fp != (NULL)) {
      char fileName[509] = "pitch.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(problem3_M, "Error closing MAT-file pitch.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(problem3_M, "Error reopening MAT-file pitch.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 2, problem3_DW.ToFile_IWORK.Count,
           "measured_pitch")) {
        rtmSetErrorStatus(problem3_M,
                          "Error writing header for measured_pitch to MAT-file pitch.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(problem3_M, "Error closing MAT-file pitch.mat");
        return;
      }

      problem3_DW.ToFile_PWORK.FilePtr = (NULL);
    }
  }

  /* Terminate for ToFile: '<Root>/To File1' */
  {
    FILE *fp = (FILE *) problem3_DW.ToFile1_PWORK.FilePtr;
    if (fp != (NULL)) {
      char fileName[509] = "travel.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(problem3_M, "Error closing MAT-file travel.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(problem3_M, "Error reopening MAT-file travel.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 2, problem3_DW.ToFile1_IWORK.Count,
           "measured_travel")) {
        rtmSetErrorStatus(problem3_M,
                          "Error writing header for measured_travel to MAT-file travel.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(problem3_M, "Error closing MAT-file travel.mat");
        return;
      }

      problem3_DW.ToFile1_PWORK.FilePtr = (NULL);
    }
  }

  /* Terminate for ToFile: '<Root>/To File2' */
  {
    FILE *fp = (FILE *) problem3_DW.ToFile2_PWORK.FilePtr;
    if (fp != (NULL)) {
      char fileName[509] = "pitch_ref.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(problem3_M, "Error closing MAT-file pitch_ref.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(problem3_M, "Error reopening MAT-file pitch_ref.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 2, problem3_DW.ToFile2_IWORK.Count,
           "pitch_ref")) {
        rtmSetErrorStatus(problem3_M,
                          "Error writing header for pitch_ref to MAT-file pitch_ref.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(problem3_M, "Error closing MAT-file pitch_ref.mat");
        return;
      }

      problem3_DW.ToFile2_PWORK.FilePtr = (NULL);
    }
  }

  /* Terminate for ToFile: '<Root>/To File3' */
  {
    FILE *fp = (FILE *) problem3_DW.ToFile3_PWORK.FilePtr;
    if (fp != (NULL)) {
      char fileName[509] = "travel_ref.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(problem3_M, "Error closing MAT-file travel_ref.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(problem3_M, "Error reopening MAT-file travel_ref.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 2, problem3_DW.ToFile3_IWORK.Count,
           "travel_ref")) {
        rtmSetErrorStatus(problem3_M,
                          "Error writing header for travel_ref to MAT-file travel_ref.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(problem3_M, "Error closing MAT-file travel_ref.mat");
        return;
      }

      problem3_DW.ToFile3_PWORK.FilePtr = (NULL);
    }
  }
}

/*========================================================================*
 * Start of Classic call interface                                        *
 *========================================================================*/

/* Solver interface called by GRT_Main */
#ifndef USE_GENERATED_SOLVER

void rt_ODECreateIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEDestroyIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEUpdateContinuousStates(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

#endif

void MdlOutputs(int_T tid)
{
  problem3_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  problem3_update();
  UNUSED_PARAMETER(tid);
}

void MdlInitializeSizes(void)
{
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
}

void MdlStart(void)
{
  problem3_initialize();
}

void MdlTerminate(void)
{
  problem3_terminate();
}

/* Registration function */
RT_MODEL_problem3_T *problem3(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  problem3_P.Integrator_UpperSat = rtInf;
  problem3_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)problem3_M, 0,
                sizeof(RT_MODEL_problem3_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&problem3_M->solverInfo,
                          &problem3_M->Timing.simTimeStep);
    rtsiSetTPtr(&problem3_M->solverInfo, &rtmGetTPtr(problem3_M));
    rtsiSetStepSizePtr(&problem3_M->solverInfo, &problem3_M->Timing.stepSize0);
    rtsiSetdXPtr(&problem3_M->solverInfo, &problem3_M->ModelData.derivs);
    rtsiSetContStatesPtr(&problem3_M->solverInfo, (real_T **)
                         &problem3_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&problem3_M->solverInfo,
      &problem3_M->Sizes.numContStates);
    rtsiSetErrorStatusPtr(&problem3_M->solverInfo, (&rtmGetErrorStatus
      (problem3_M)));
    rtsiSetRTModelPtr(&problem3_M->solverInfo, problem3_M);
  }

  rtsiSetSimTimeStep(&problem3_M->solverInfo, MAJOR_TIME_STEP);
  problem3_M->ModelData.intgData.f[0] = problem3_M->ModelData.odeF[0];
  problem3_M->ModelData.contStates = ((real_T *) &problem3_X);
  rtsiSetSolverData(&problem3_M->solverInfo, (void *)
                    &problem3_M->ModelData.intgData);
  rtsiSetSolverName(&problem3_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = problem3_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    problem3_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    problem3_M->Timing.sampleTimes = (&problem3_M->Timing.sampleTimesArray[0]);
    problem3_M->Timing.offsetTimes = (&problem3_M->Timing.offsetTimesArray[0]);

    /* task periods */
    problem3_M->Timing.sampleTimes[0] = (0.0);
    problem3_M->Timing.sampleTimes[1] = (0.002);

    /* task offsets */
    problem3_M->Timing.offsetTimes[0] = (0.0);
    problem3_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(problem3_M, &problem3_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = problem3_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    problem3_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(problem3_M, -1);
  problem3_M->Timing.stepSize0 = 0.002;
  problem3_M->Timing.stepSize1 = 0.002;

  /* External mode info */
  problem3_M->Sizes.checksums[0] = (1486277045U);
  problem3_M->Sizes.checksums[1] = (1408228525U);
  problem3_M->Sizes.checksums[2] = (1591682926U);
  problem3_M->Sizes.checksums[3] = (191186349U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    problem3_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(problem3_M->extModeInfo,
      &problem3_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(problem3_M->extModeInfo, problem3_M->Sizes.checksums);
    rteiSetTPtr(problem3_M->extModeInfo, rtmGetTPtr(problem3_M));
  }

  problem3_M->solverInfoPtr = (&problem3_M->solverInfo);
  problem3_M->Timing.stepSize = (0.002);
  rtsiSetFixedStepSize(&problem3_M->solverInfo, 0.002);
  rtsiSetSolverMode(&problem3_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  problem3_M->ModelData.blockIO = ((void *) &problem3_B);

  {
    problem3_B.TravelCounttorad = 0.0;
    problem3_B.Gain = 0.0;
    problem3_B.Sum3 = 0.0;
    problem3_B.Gain1 = 0.0;
    problem3_B.PitchCounttorad = 0.0;
    problem3_B.Gain_i = 0.0;
    problem3_B.Gain1_f = 0.0;
    problem3_B.Gain_d = 0.0;
    problem3_B.Gain_b = 0.0;
    problem3_B.Sum4 = 0.0;
    problem3_B.ElevationCounttorad = 0.0;
    problem3_B.Gain_e = 0.0;
    problem3_B.Sum = 0.0;
    problem3_B.Gain_dg = 0.0;
    problem3_B.Sum1 = 0.0;
    problem3_B.Sum2 = 0.0;
    problem3_B.K_ei = 0.0;
    problem3_B.Gain_l = 0.0;
    problem3_B.BackmotorSaturation = 0.0;
    problem3_B.FrontmotorSaturation = 0.0;
    problem3_B.FromWorkspace = 0.0;
  }

  /* parameters */
  problem3_M->ModelData.defaultParam = ((real_T *)&problem3_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &problem3_X;
    problem3_M->ModelData.contStates = (x);
    (void) memset((void *)&problem3_X, 0,
                  sizeof(X_problem3_T));
  }

  /* states (dwork) */
  problem3_M->ModelData.dwork = ((void *) &problem3_DW);
  (void) memset((void *)&problem3_DW, 0,
                sizeof(DW_problem3_T));

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      problem3_DW.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      problem3_DW.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      problem3_DW.HILInitialize_AOMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      problem3_DW.HILInitialize_AOMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      problem3_DW.HILInitialize_AOVoltages[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      problem3_DW.HILInitialize_FilterFrequency[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      problem3_DW.HILInitialize_POSortedFreqs[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      problem3_DW.HILInitialize_POValues[i] = 0.0;
    }
  }

  problem3_DW.TimeStampA = 0.0;
  problem3_DW.LastUAtTimeA = 0.0;
  problem3_DW.TimeStampB = 0.0;
  problem3_DW.LastUAtTimeB = 0.0;
  problem3_DW.HILWriteAnalog_Buffer[0] = 0.0;
  problem3_DW.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    problem3_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 16;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.B = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.P = &rtPTransTable;
  }

  /* Initialize Sizes */
  problem3_M->Sizes.numContStates = (4);/* Number of continuous states */
  problem3_M->Sizes.numY = (0);        /* Number of model outputs */
  problem3_M->Sizes.numU = (0);        /* Number of model inputs */
  problem3_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  problem3_M->Sizes.numSampTimes = (2);/* Number of sample times */
  problem3_M->Sizes.numBlocks = (69);  /* Number of blocks */
  problem3_M->Sizes.numBlockIO = (21); /* Number of block outputs */
  problem3_M->Sizes.numBlockPrms = (149);/* Sum of parameter "widths" */
  return problem3_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
