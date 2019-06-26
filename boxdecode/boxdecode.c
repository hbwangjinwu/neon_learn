#include"sample_svp_nnie_software.h"
#include <math.h>
#include <sys/time.h>
#include  <stdio.h>
#include <omp.h>

#define  threads_num 4

#ifdef __cplusplus    // If used by C++ code,
extern "C" {          // we need to export the C interface
#endif
static struct timeval t_start, t_end;
static int timediff(struct timeval a, struct timeval b)
{
    return (a.tv_sec - b.tv_sec) * 1000000 + (a.tv_usec - b.tv_usec);
}
#define time_start() \
    gettimeofday(&t_start, 0);
	
#define time_stop(tag)          \
    gettimeofday(&t_end, 0); \
    printf("%s : %5d us\n", tag, timediff(t_end, t_start));

#define time_reset(tag) \
    time_stop(tag);     \
    gettimeofday(&t_start, 0);

void calDetectBox(HI_S32* prior,HI_S32* Priorvar,HI_S32* locpred,HI_S32* decodebox )
{
	
	HI_FLOAT f32PriorWidth = (HI_FLOAT)(prior[2] - prior[0]);
	HI_FLOAT f32PriorHeight = (HI_FLOAT)(prior[3] - prior[1]);
	HI_FLOAT f32PriorCenterX = (prior[2] + prior[0])*SAMPLE_SVP_NNIE_HALF;
	HI_FLOAT f32PriorCenterY = (prior[3] + prior[1])*SAMPLE_SVP_NNIE_HALF;

	HI_FLOAT f32DecodeBoxCenterX = ((HI_FLOAT)Priorvar[0]/SAMPLE_SVP_NNIE_QUANT_BASE)*
		((HI_FLOAT)locpred[0]/SAMPLE_SVP_NNIE_QUANT_BASE)*f32PriorWidth+f32PriorCenterX;

	HI_FLOAT f32DecodeBoxCenterY = ((HI_FLOAT)Priorvar[1]/SAMPLE_SVP_NNIE_QUANT_BASE)*
		((HI_FLOAT)locpred[1]/SAMPLE_SVP_NNIE_QUANT_BASE)*f32PriorHeight+f32PriorCenterY;
		
	HI_FLOAT f32DecodeBoxWidth = SVP_NNIE_QuickExp((HI_S32)(((HI_FLOAT)Priorvar[2]/SAMPLE_SVP_NNIE_QUANT_BASE)*
		locpred[2]))*f32PriorWidth;
	HI_FLOAT f32DecodeBoxHeight = SVP_NNIE_QuickExp((HI_S32)(((HI_FLOAT)Priorvar[3]/SAMPLE_SVP_NNIE_QUANT_BASE)*
		locpred[3]))*f32PriorHeight;
	decodebox[0] = (HI_S32)(f32DecodeBoxCenterX - f32DecodeBoxWidth * SAMPLE_SVP_NNIE_HALF);
	decodebox[1] = (HI_S32)(f32DecodeBoxCenterY - f32DecodeBoxHeight * SAMPLE_SVP_NNIE_HALF);
	decodebox[2] = (HI_S32)(f32DecodeBoxCenterX + f32DecodeBoxWidth * SAMPLE_SVP_NNIE_HALF);
	decodebox[3] = (HI_S32)(f32DecodeBoxCenterY + f32DecodeBoxHeight * SAMPLE_SVP_NNIE_HALF);
}


void calDetectBoxNEON(HI_S32* prior,HI_S32* Priorvar,HI_S32* locpred,HI_S32* decodebox )
{
	float half=0.5;
	float quant=4096.0;
	asm volatile(
		"dup    v0.4s,%w[half]                     \n"
		"dup    v1.4s,%w[quant_base]               \n"
		"frecpe v3.4s,v1.4s                        \n"
		"frecps v2.4s,v1.4s,v3.4s                  \n"
		"fmul   v1.4s,v2.4s,v3.4s                  \n" 
		
		"prfm   pldl1keep, [%[prioraddr], #128] 		\n"
		"ld4    {v2.4s,v3.4s,v4.4s,v5.4s},[%[prioraddr]] \n"
		"scvtf  v2.4s,v2.4s                      \n"
		"scvtf  v3.4s,v3.4s                      \n"
		"scvtf  v4.4s,v4.4s                      \n"
		"scvtf  v5.4s,v5.4s                      \n"
		"fsub   v6.4s,v4.4s,v2.4s   			 \n" //f32PriorWidth
		"fsub 	v7.4s,v5.4s,v3.4s                \n" //f32PriorHeight
		"fadd	v8.4s,v4.4s,v2.4s				 \n"
		"fadd 	v9.4s,v5.4s,v3.4s				 \n"
		"fmul	v8.4s,v8.4s,v0.4s                \n" //f32PriorCenterX
		"fmul 	v9.4s,v9.4s,v0.4s				 \n" //f32PriorCenterY              .0
		
		"prfm   pldl1keep, [%[priorvaraddr], #128] 		\n"
		"ld4    {v10.4s,v11.4s,v12.4s,v13.4s},[%[priorvaraddr]] \n"
		"scvtf  v10.4s,v10.4s                      \n" 
		"scvtf  v11.4s,v11.4s                      \n" 
		"scvtf  v12.4s,v12.4s                      \n" 
		"scvtf  v13.4s,v13.4s                      \n" 
		
		"fmul   v10.4s,v10.4s,v1.4s                      \n"//Priorvar[0]/SAMPLE_SVP_NNIE_QUANT_BASE
		"fmul   v11.4s,v11.4s,v1.4s                      \n"//Priorvar[1]/SAMPLE_SVP_NNIE_QUANT_BASE
		"fmul   v12.4s,v12.4s,v1.4s                      \n"//Priorvar[2]/SAMPLE_SVP_NNIE_QUANT_BASE
		"fmul   v13.4s,v13.4s,v1.4s                      \n"//Priorvar[3]/SAMPLE_SVP_NNIE_QUANT_BASE 
		
		
		"prfm   pldl1keep, [%[locpredaddr], #128] 		\n"
		"ld4    {v14.4s,v15.4s,v16.4s,v17.4s},[%[locpredaddr]] \n"
		"scvtf  v14.4s,v14.4s                      \n"
		"scvtf  v15.4s,v15.4s                      \n"
		"scvtf  v16.4s,v16.4s                      \n"
		"scvtf  v17.4s,v17.4s                      \n"
		
		"fmul   v14.4s,v14.4s,v1.4s                      \n"//locpred[0]/SAMPLE_SVP_NNIE_QUANT_BASE
		"fmul   v15.4s,v15.4s,v1.4s                      \n"//ocpred[1]/SAMPLE_SVP_NNIE_QUANT_BASE
		
		"fmul   v2.4s,v10.4s,v14.4s                  \n"
		"fmul   v3.4s,v11.4s,v15.4s                  \n"
		"fmul   v2.4s,v2.4s,v6.4s                    \n"
		"fmul   v3.4s,v3.4s,v7.4s                    \n"
		"fadd   v2.4s,v2.4s,v8.4s                    \n" //f32DecodeBoxCenterX
		"fadd   v3.4s,v3.4s,v9.4s                    \n" //f32DecodeBoxCenterY
		
		//v10 v11 v14  v15  now useable 
		"fmul   v4.4s,v12.4s,v16.4s                  \n"
		"fmul   v5.4s,v13.4s,v17.4s					 \n"
		"fcvtzs v4.4s,v4.4s                          \n"  //(HI_S32)(((HI_FLOAT)Priorvar[2]/SAMPLE_SVP_NNIE_QUANT_BASE)*locpred[2])
		"fcvtzs v5.4s,v5.4s                          \n"  //HI_S32)(((HI_FLOAT)Priorvar[3]/SAMPLE_SVP_NNIE_QUANT_BASE)*locpred[3])
		//v4 v5  value is positive ----need be sure 
		
		
		
		//只处理了正数 待修正
		/*
		"mov  	w0,#0x3ff                             \n"
		"dup    v10.4s,w0                             \n"
		"and    v14.16b,v4.16b,v10.16b                \n" //v4&0x3ff
		"and    v15.16b,v5.16b,v10.16b                \n" //v5&0x3ff
		
		"sshr   v4.4s,v4.4s,#10                       \n"
		"sshr   v5.4s,v5.4s,#10                       \n"
		"and    v4.16b,v4.16b,v10.16b                \n" //(v4>>10)&0x3ff
		"and    v5.16b,v5.16b,v10.16b                \n" //(v5>>10)&0x3ff
		"mov    w0,#1024                             \n"
		"dup    v10.4s,w0                            \n"
		"add    v4.4s,v4.4s,v10.4s                   \n"
		"add    v5.4s,v5.4s,v10.4s                   \n"
		*/
		//以下为修正代码
		"mov  	w0,#0                                 \n"
		"dup  	v10.4s,w0                                \n"
		"cmgt   v11.4s,v4.4s,v10.4s\n"   //v11 select 
		"cmgt   v12.4s,v5.4s,v10.4s\n"   //v12 select
		
		"abs    v4.4s,v4.4s    \n"   
		"abs    v5.4s,v5.4s    \n"		
		//low 
		"mov    v14.16b,v11.16b   \n"
		"mov    v15.16b,v12.16b   \n"
		//"mov    w0,#0           \n"
		//"dup    v10.4s,w0         \n"
		"mov    w0,#2048        \n"
		"dup    v18.4s,w0         \n"
		"bsl    v14.16b,v10.16b,v18.16b  \n" //  low	
		"bsl    v15.16b,v10.16b,v18.16b  \n" //  low	
		
		"mov    w0,#0x3ff        \n"
		"dup    v10.4s,w0         \n"
		"and    v16.16b, v4.16b, v10.16b \n" //
		"and    v17.16b, v5.16b, v10.16b \n" //
		"add    v14.4s,v14.4s, v16.4s    \n"
		"add    v15.4s,v15.4s, v17.4s    \n"
		//high
		"mov    v18.16b,v11.16b   \n"
		"mov    v19.16b,v12.16b   \n"
		
		"mov    w0,#1024           \n"
		"dup    v16.4s,w0         \n"
		"mov    w0,#3072        \n"
		"dup    v17.4s,w0         \n"
		
		"bsl    v18.16b,v16.16b,v17.16b  \n" //
		"bsl    v19.16b,v16.16b,v17.16b  \n"
		"sshr   v4.4s, v4.4s,#10 \n"
		"sshr   v5.4s, v5.4s,#10 \n"
		"mov    w0,#0x3ff        \n"
		"dup    v16.4s,w0         \n"
		"and    v10.16b, v4.16b, v16.16b \n" //
		"and    v11.16b, v5.16b, v16.16b \n" //
		"add    v4.4s,v18.4s, v10.4s    \n"
		"add    v5.4s,v19.4s, v11.4s    \n"
		
		"mov    x0,#0                                \n" //make sure the whole offset for coef is right[low is reassigned by w0]
		"mov    w0,v14.s[0]                          \n"
		"ldr    w0,[%[expcoef],x0,lsl#2]	         \n"
		"mov    v10.s[0],w0                          \n"
		
		"mov    w0,v15.s[0]                          \n"
		"ldr    w0,[%[expcoef],x0,lsl#2]	         \n"
		"mov    v12.s[0],w0                          \n"
		
		"mov    w0,v14.s[1]                          \n"
		"ldr    w0,[%[expcoef],x0,lsl#2]	         \n"
		"mov    v10.s[1],w0                          \n"
		
		"mov    w0,v15.s[1]                          \n"
		"ldr    w0,[%[expcoef],x0,lsl#2]	         \n"
		"mov    v12.s[1],w0                          \n"
		
		"mov    w0,v14.s[2]                          \n"
		"ldr    w0,[%[expcoef],x0,lsl#2]	         \n"
		"mov    v10.s[2],w0                          \n"
		
		"mov    w0,v15.s[2]                          \n"
		"ldr    w0,[%[expcoef],x0,lsl#2]	         \n"
		"mov    v12.s[2],w0                          \n"
		
		"mov    w0,v14.s[3]                          \n"
		"ldr    w0,[%[expcoef],x0,lsl#2]	         \n"
		"mov    v10.s[3],w0                          \n"
		
		"mov    w0,v15.s[3]                          \n"
		"ldr    w0,[%[expcoef],x0,lsl#2]	         \n"
		"mov    v12.s[3],w0                          \n"
		
		"mov    w0,v4.s[0]                          \n"
		"ldr    w0,[%[expcoef],x0,lsl#2]	         \n"
		"mov    v11.s[0],w0                          \n"
		
		"mov    w0,v5.s[0]                           \n"
		"ldr    w0,[%[expcoef],x0,lsl#2]	         \n"
		"mov    v13.s[0],w0                          \n"
		
		"mov    w0,v4.s[1]                          \n"
		"ldr    w0,[%[expcoef],x0,lsl#2]	         \n"
		"mov    v11.s[1],w0                          \n"
		
		"mov    w0,v5.s[1]                          \n"
		"ldr    w0,[%[expcoef],x0,lsl#2]	         \n"
		"mov    v13.s[1],w0                          \n"
		
		"mov    w0,v4.s[2]                          \n"
		"ldr    w0,[%[expcoef],x0,lsl#2]	         \n"
		"mov    v11.s[2],w0                          \n"
		
		"mov    w0,v5.s[2]                          \n"
		"ldr    w0,[%[expcoef],x0,lsl#2]	         \n"
		"mov    v13.s[2],w0                          \n"
		
		"mov    w0,v4.s[3]                          \n"
		"ldr    w0,[%[expcoef],x0,lsl#2]	         \n"
		"mov    v11.s[3],w0                          \n"
		
		"mov    w0,v5.s[3]                          \n"
		"ldr    w0,[%[expcoef],x0,lsl#2]	         \n"
		"mov    v13.s[3],w0                          \n"
		
		"fmul  v4.4s,v10.4s,v11.4s                   \n"//SVP_NNIE_QuickExp((HI_S32)(((HI_FLOAT)Priorvar[2]/SAMPLE_SVP_NNIE_QUANT_BASE)*locpred[2]))
		"fmul  v5.4s,v12.4s,v13.4s                   \n"//SVP_NNIE_QuickExp((HI_S32)(((HI_FLOAT)Priorvar[3]/SAMPLE_SVP_NNIE_QUANT_BASE)*locpred[3]))
		
		"fmul  v4.4s,v4.4s,v6.4s                     \n"//f32DecodeBoxWidth
		"fmul  v5.4s,v5.4s,v7.4s				     \n"//f32DecodeBoxHeight
		
		"fmul v10.4s,v4.4s,v0.4s   \n" //f32DecodeBoxWidth * SAMPLE_SVP_NNIE_HALF
		"fmul v11.4s,v5.4s,v0.4s   \n" //f32DecodeBoxHeight * SAMPLE_SVP_NNIE_HALF
		
		"fsub v6.4s,v2.4s,v10.4s   \n"
		"fsub v7.4s,v3.4s,v11.4s   \n"
		"fadd v8.4s,v2.4s,v10.4s   \n"
		"fadd v9.4s,v3.4s,v11.4s   \n"
		
		"fcvtzs v6.4s,v6.4s        \n"
		"fcvtzs v7.4s,v7.4s        \n"
		"fcvtzs v8.4s,v8.4s        \n"
		"fcvtzs v9.4s,v9.4s        \n"
		"st4  {v6.4s,v7.4s,v8.4s,v9.4s},[%[decodeboxaddr]] \n" 
		:
		:[prioraddr]"r"(prior),
		 [priorvaraddr]"r"(Priorvar),
		 [locpredaddr]"r"(locpred),
		 [decodeboxaddr]"r"(decodebox),
		 [half]"r"(half),
		 [expcoef]"r"(s_af32ExpCoef),
		 [quant_base]"r"(quant)
		:"cc","memory","w0","v0","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","v11","v12","v13","v14","v15","v16","v17","v18","v19"
	
	);
}


/*****************************************************************************
* Prototype :   SVP_NNIE_Ssd_DetectionOutForward
* Description : this function is used to get detection result of SSD
* Input :     HI_U32 u32ConcatNum            [IN] SSD concat num
*             HI_U32 u32ConfThresh           [IN] confidence thresh
*             HI_U32 u32ClassNum             [IN] class num
*             HI_U32 u32TopK                 [IN] Topk value
*             HI_U32 u32KeepTopK             [IN] KeepTopK value
*             HI_U32 u32NmsThresh            [IN] NMS thresh
*             HI_U32 au32DetectInputChn[]    [IN] detection input channel
*             HI_S32* aps32AllLocPreds[]     [IN] Location prediction
*             HI_S32* aps32AllPriorBoxes[]   [IN] prior box
*             HI_S32* ps32ConfScores         [IN] confidence score
*             HI_S32* ps32AssistMemPool      [IN] assist buffer
*             HI_S32* ps32DstScoreSrc        [OUT] result of score
*             HI_S32* ps32DstBboxSrc         [OUT] result of Bbox
*             HI_S32* ps32RoiOutCntSrc       [OUT] result of the roi num of each class
*
*
* Output :
* Return Value : HI_SUCCESS: Success;Error codes: Failure.
* Spec :
* Calls :
* Called By :
* History:
*
* 1. Date : 2017-11-10
* Author :
* Modification : Create
*
*****************************************************************************/
static HI_S32 SVP_NNIE_Ssd_DetectionOutForward(HI_U32 u32ConcatNum,
    HI_U32 u32ConfThresh,HI_U32 u32ClassNum, HI_U32 u32TopK, HI_U32 u32KeepTopK, HI_U32 u32NmsThresh,
    HI_U32 au32DetectInputChn[], HI_S32* aps32AllLocPreds[], HI_S32* aps32AllPriorBoxes[],
    HI_S32* ps32ConfScores, HI_S32* ps32AssistMemPool, HI_S32* ps32DstScoreSrc,
    HI_S32* ps32DstBboxSrc, HI_S32* ps32RoiOutCntSrc)
{
    HI_S32* ps32LocPreds = NULL;
    HI_S32* ps32PriorBoxes = NULL;
    HI_S32* ps32PriorVar = NULL;
    HI_S32* ps32AllDecodeBoxes = NULL;
    HI_S32* ps32DstScore = NULL;
    HI_S32* ps32DstBbox = NULL;
    HI_S32* ps32ClassRoiNum = NULL;
    HI_U32 u32RoiOutCnt = 0;
    HI_S32* ps32SingleProposal = NULL;
    HI_S32* ps32AfterTopK = NULL;
    SAMPLE_SVP_NNIE_STACK_S* pstStack = NULL;
    HI_U32 u32PriorNum = 0;
    HI_U32 u32NumPredsPerClass = 0;
    HI_U32 u32SrcIdx = 0;
    HI_U32 u32AfterFilter = 0;
    HI_U32 u32AfterTopK = 0;
    HI_U32 u32AfterConf = 0;
    HI_U32 u32KeepCnt = 0;
    HI_U32 i = 0;
    HI_U32 j = 0;
    HI_U32 u32Offset = 0;
    HI_S32 s32Ret = HI_SUCCESS;
    u32PriorNum = 0;
	HI_U32 index=0;
    for (i = 0; i < u32ConcatNum; i++)
    {
        u32PriorNum += au32DetectInputChn[i] / SAMPLE_SVP_NNIE_COORDI_NUM;
    }
    ps32AllDecodeBoxes = ps32AssistMemPool;
    ps32SingleProposal = ps32AllDecodeBoxes + u32PriorNum * SAMPLE_SVP_NNIE_COORDI_NUM;
    ps32AfterTopK = ps32SingleProposal + SAMPLE_SVP_NNIE_PROPOSAL_WIDTH * u32PriorNum;
    pstStack = (SAMPLE_SVP_NNIE_STACK_S*)(ps32AfterTopK + u32PriorNum * SAMPLE_SVP_NNIE_PROPOSAL_WIDTH);
	//time_start();
#if 0
	HI_S32* DecodeBox=ps32AllDecodeBoxes;
	u32SrcIdx = 0;
    for (i = 0; i < u32ConcatNum; i++)
    {
        /********** get loc predictions ************/
        ps32LocPreds = aps32AllLocPreds[i];
        u32NumPredsPerClass = au32DetectInputChn[i] / SAMPLE_SVP_NNIE_COORDI_NUM;
        /********** get Prior Bboxes ************/
        ps32PriorBoxes = aps32AllPriorBoxes[i];
        ps32PriorVar = ps32PriorBoxes + u32NumPredsPerClass*SAMPLE_SVP_NNIE_COORDI_NUM;
		#pragma omp parallel for num_threads(threads_num)
		for (j = 0; j < u32NumPredsPerClass; j++)
        {		
			HI_S32* prior=ps32PriorBoxes+j*SAMPLE_SVP_NNIE_COORDI_NUM;
			HI_S32* Priorvar=ps32PriorVar+j*SAMPLE_SVP_NNIE_COORDI_NUM;
			HI_S32* locpred=ps32LocPreds+j*SAMPLE_SVP_NNIE_COORDI_NUM;
			HI_S32* decodebox=DecodeBox+j*SAMPLE_SVP_NNIE_COORDI_NUM;
			calDetectBox(prior,Priorvar,locpred,decodebox );
		}
		DecodeBox+=SAMPLE_SVP_NNIE_COORDI_NUM*u32NumPredsPerClass;
    }
#else 
	HI_S32* DecodeBox=ps32AllDecodeBoxes;
	u32SrcIdx = 0;
    for (i = 0; i < u32ConcatNum; i++)
    {
        /********** get loc predictions ************/
        ps32LocPreds = aps32AllLocPreds[i];
        u32NumPredsPerClass = au32DetectInputChn[i] / SAMPLE_SVP_NNIE_COORDI_NUM;
        /********** get Prior Bboxes ************/
        ps32PriorBoxes = aps32AllPriorBoxes[i];
        ps32PriorVar = ps32PriorBoxes + u32NumPredsPerClass*SAMPLE_SVP_NNIE_COORDI_NUM;
		HI_S32 nn=u32NumPredsPerClass>>2;
		HI_S32 remain_start=nn<<2;
		HI_S32 remain=u32NumPredsPerClass-remain_start;
		
		#pragma omp parallel for num_threads(threads_num)
		for(j=0;j<nn;j++)
		{
			HI_S32* prior=ps32PriorBoxes+j*SAMPLE_SVP_NNIE_COORDI_NUM*4;
			HI_S32* Priorvar=ps32PriorVar+j*SAMPLE_SVP_NNIE_COORDI_NUM*4;
			HI_S32* locpred=ps32LocPreds+j*SAMPLE_SVP_NNIE_COORDI_NUM*4;
			HI_S32* decodebox=DecodeBox+j*SAMPLE_SVP_NNIE_COORDI_NUM*4;
			calDetectBoxNEON(prior,Priorvar,locpred,decodebox );
		}
		for (j = remain_start; j < u32NumPredsPerClass; j++)
        {		
			HI_S32* prior=ps32PriorBoxes+j*SAMPLE_SVP_NNIE_COORDI_NUM;
			HI_S32* Priorvar=ps32PriorVar+j*SAMPLE_SVP_NNIE_COORDI_NUM;
			HI_S32* locpred=ps32LocPreds+j*SAMPLE_SVP_NNIE_COORDI_NUM;
			HI_S32* decodebox=DecodeBox+j*SAMPLE_SVP_NNIE_COORDI_NUM;
			calDetectBox(prior,Priorvar,locpred,decodebox );
		}
		DecodeBox+=SAMPLE_SVP_NNIE_COORDI_NUM*u32NumPredsPerClass;
    }
#endif
	//time_stop("detect box ");
    /********** do NMS for each class *************/
    HI_U32 u32NonBackground = 0;
	u32ConfThresh=1024;
    for (j = 0; j < u32PriorNum; j++)
    {
		if(ps32ConfScores[j*u32ClassNum + 0]>u32ConfThresh)continue;
		ps32AfterTopK[u32NonBackground * SAMPLE_SVP_NNIE_PROPOSAL_WIDTH + 0] = ps32AllDecodeBoxes[j * SAMPLE_SVP_NNIE_COORDI_NUM + 0];
		ps32AfterTopK[u32NonBackground * SAMPLE_SVP_NNIE_PROPOSAL_WIDTH + 1] = ps32AllDecodeBoxes[j * SAMPLE_SVP_NNIE_COORDI_NUM + 1];
		ps32AfterTopK[u32NonBackground * SAMPLE_SVP_NNIE_PROPOSAL_WIDTH + 2] = ps32AllDecodeBoxes[j * SAMPLE_SVP_NNIE_COORDI_NUM + 2];
		ps32AfterTopK[u32NonBackground * SAMPLE_SVP_NNIE_PROPOSAL_WIDTH + 3] = ps32AllDecodeBoxes[j * SAMPLE_SVP_NNIE_COORDI_NUM + 3];
		//ps32AfterTopK[u32NonBackground * SAMPLE_SVP_NNIE_PROPOSAL_WIDTH + 4] = ps32ConfScores[j*u32ClassNum + i];
		ps32AfterTopK[u32NonBackground * SAMPLE_SVP_NNIE_PROPOSAL_WIDTH + 5] = j;
		u32NonBackground++;
    }
	u32PriorNum = u32NonBackground;
	//printf("u32NonBG = %d\n",u32NonBackground);
	u32AfterTopK = 0;
	for (i = 1; i < u32ClassNum; i++)
	{
		u32AfterConf = 0;
		for (j = 0; j < u32PriorNum; j++)
		{
			HI_S32 s32Score = ps32ConfScores[ps32AfterTopK[j * SAMPLE_SVP_NNIE_PROPOSAL_WIDTH + 5]*u32ClassNum + i];
			if(s32Score > u32NmsThresh){
				ps32SingleProposal[u32AfterConf * SAMPLE_SVP_NNIE_PROPOSAL_WIDTH + 0] = ps32AfterTopK[j*SAMPLE_SVP_NNIE_PROPOSAL_WIDTH + 0];
				ps32SingleProposal[u32AfterConf * SAMPLE_SVP_NNIE_PROPOSAL_WIDTH + 1] = ps32AfterTopK[j*SAMPLE_SVP_NNIE_PROPOSAL_WIDTH + 1];
				ps32SingleProposal[u32AfterConf * SAMPLE_SVP_NNIE_PROPOSAL_WIDTH + 2] = ps32AfterTopK[j*SAMPLE_SVP_NNIE_PROPOSAL_WIDTH + 2];
				ps32SingleProposal[u32AfterConf * SAMPLE_SVP_NNIE_PROPOSAL_WIDTH + 3] = ps32AfterTopK[j*SAMPLE_SVP_NNIE_PROPOSAL_WIDTH + 3];
				ps32SingleProposal[u32AfterConf * SAMPLE_SVP_NNIE_PROPOSAL_WIDTH + 4] = s32Score;
				ps32SingleProposal[u32AfterConf * SAMPLE_SVP_NNIE_PROPOSAL_WIDTH + 5] = 0;
				u32AfterConf ++;
			}
        }
        s32Ret = SVP_NNIE_NonRecursiveArgQuickSort(ps32SingleProposal, 0, u32AfterConf - 1, pstStack, u32TopK);
        u32AfterFilter = (u32AfterConf < u32TopK) ? u32AfterConf : u32TopK;
        s32Ret = SVP_NNIE_NonMaxSuppression(ps32SingleProposal, u32AfterFilter, u32NmsThresh, u32AfterFilter);
        u32RoiOutCnt = 0;
        ps32DstScore = (HI_S32*)ps32DstScoreSrc;
        ps32DstBbox = (HI_S32*)ps32DstBboxSrc;
        ps32ClassRoiNum = (HI_S32*)ps32RoiOutCntSrc;
        ps32DstScore += (HI_S32)u32AfterTopK;
        ps32DstBbox += (HI_S32)(u32AfterTopK * SAMPLE_SVP_NNIE_COORDI_NUM);
        for (j = 0; j < u32AfterFilter; j++)
        {
            if (ps32SingleProposal[j * SAMPLE_SVP_NNIE_PROPOSAL_WIDTH + 5] == 0)
            {
                ps32DstScore[u32RoiOutCnt] = ps32SingleProposal[j * 6 + 4];
                ps32DstBbox[u32RoiOutCnt * SAMPLE_SVP_NNIE_COORDI_NUM] = ps32SingleProposal[j * SAMPLE_SVP_NNIE_PROPOSAL_WIDTH];
                ps32DstBbox[u32RoiOutCnt * SAMPLE_SVP_NNIE_COORDI_NUM + 1] = ps32SingleProposal[j * SAMPLE_SVP_NNIE_PROPOSAL_WIDTH + 1];
                ps32DstBbox[u32RoiOutCnt * SAMPLE_SVP_NNIE_COORDI_NUM + 2] = ps32SingleProposal[j * SAMPLE_SVP_NNIE_PROPOSAL_WIDTH + 2];
                ps32DstBbox[u32RoiOutCnt * SAMPLE_SVP_NNIE_COORDI_NUM + 3] = ps32SingleProposal[j * SAMPLE_SVP_NNIE_PROPOSAL_WIDTH + 3];
                u32RoiOutCnt++;
            }
        }
        ps32ClassRoiNum[i] = (HI_S32)u32RoiOutCnt;
        u32AfterTopK += u32RoiOutCnt;
	    // time_stop("---->post");
    }
    //time_reset("-->nms");
    u32KeepCnt = 0;
    u32Offset = 0;
    if (u32AfterTopK > u32KeepTopK)
    {
        u32Offset = ps32ClassRoiNum[0];
        for (i = 1; i < u32ClassNum; i++)
        {
            ps32DstScore = (HI_S32*)ps32DstScoreSrc;
            ps32DstBbox = (HI_S32*)ps32DstBboxSrc;
            ps32ClassRoiNum = (HI_S32*)ps32RoiOutCntSrc;
            ps32DstScore += (HI_S32)(u32Offset);
            ps32DstBbox += (HI_S32)(u32Offset * SAMPLE_SVP_NNIE_COORDI_NUM);
            for (j = 0; j < (HI_U32)ps32ClassRoiNum[i]; j++)
            {
                ps32AfterTopK[u32KeepCnt * SAMPLE_SVP_NNIE_PROPOSAL_WIDTH] = ps32DstBbox[j * SAMPLE_SVP_NNIE_COORDI_NUM];
                ps32AfterTopK[u32KeepCnt * SAMPLE_SVP_NNIE_PROPOSAL_WIDTH + 1] = ps32DstBbox[j * SAMPLE_SVP_NNIE_COORDI_NUM + 1];
                ps32AfterTopK[u32KeepCnt * SAMPLE_SVP_NNIE_PROPOSAL_WIDTH + 2] = ps32DstBbox[j * SAMPLE_SVP_NNIE_COORDI_NUM + 2];
                ps32AfterTopK[u32KeepCnt * SAMPLE_SVP_NNIE_PROPOSAL_WIDTH + 3] = ps32DstBbox[j * SAMPLE_SVP_NNIE_COORDI_NUM + 3];
                ps32AfterTopK[u32KeepCnt * SAMPLE_SVP_NNIE_PROPOSAL_WIDTH + 4] = ps32DstScore[j];
                ps32AfterTopK[u32KeepCnt * SAMPLE_SVP_NNIE_PROPOSAL_WIDTH + 5] = i;
                u32KeepCnt++;
            }
            u32Offset = u32Offset + ps32ClassRoiNum[i];
        }
        s32Ret = SVP_NNIE_NonRecursiveArgQuickSort(ps32AfterTopK, 0, u32KeepCnt - 1, pstStack,u32KeepCnt);

        u32Offset = 0;
        u32Offset = ps32ClassRoiNum[0];
        for (i = 1; i < u32ClassNum; i++)
        {
            u32RoiOutCnt = 0;
            ps32DstScore = (HI_S32*)ps32DstScoreSrc;
            ps32DstBbox = (HI_S32*)ps32DstBboxSrc;
            ps32ClassRoiNum = (HI_S32*)ps32RoiOutCntSrc;
            ps32DstScore += (HI_S32)(u32Offset);
            ps32DstBbox += (HI_S32)(u32Offset * SAMPLE_SVP_NNIE_COORDI_NUM);
            for (j = 0; j < u32KeepTopK; j++)
            {
                if (ps32AfterTopK[j * SAMPLE_SVP_NNIE_PROPOSAL_WIDTH + 5] == i)
                {
                    ps32DstScore[u32RoiOutCnt] = ps32AfterTopK[j * SAMPLE_SVP_NNIE_PROPOSAL_WIDTH + 4];
                    ps32DstBbox[u32RoiOutCnt * SAMPLE_SVP_NNIE_COORDI_NUM] = ps32AfterTopK[j * SAMPLE_SVP_NNIE_PROPOSAL_WIDTH];
                    ps32DstBbox[u32RoiOutCnt * SAMPLE_SVP_NNIE_COORDI_NUM + 1] = ps32AfterTopK[j * SAMPLE_SVP_NNIE_PROPOSAL_WIDTH + 1];
                    ps32DstBbox[u32RoiOutCnt * SAMPLE_SVP_NNIE_COORDI_NUM + 2] = ps32AfterTopK[j * SAMPLE_SVP_NNIE_PROPOSAL_WIDTH + 2];
                    ps32DstBbox[u32RoiOutCnt * SAMPLE_SVP_NNIE_COORDI_NUM + 3] = ps32AfterTopK[j * SAMPLE_SVP_NNIE_PROPOSAL_WIDTH + 3];
                    u32RoiOutCnt++;
                }
            }
            ps32ClassRoiNum[i] = (HI_S32)u32RoiOutCnt;
            u32Offset += u32RoiOutCnt;
        }
    }
    //time_reset("-->aftertop");
    return s32Ret;
}

#ifdef __cplusplus
}
#endif
