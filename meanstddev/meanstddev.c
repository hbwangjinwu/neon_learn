#include "yuvprocess.h"
#include "sample_comm.h"
#include "sample_comm_ive.h"
#include <math.h>
#include "sys/time.h"
#include <sys/time.h>


//#define _USE_MMZ_
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
/*******************meanstddev*******************************************
 *description: calculate mean and stddev value for a input mat one channel 
               this function is optimize with neon for CV_32FC1 type.
               convert the mat to CV_32FC1 data type inhead 
 *parameter:
 * src :the data pointer to the source data
 * widht : the width infomation of the input picture
 * height: the height information of the input picture
 * mean: pointer to the mean value
 * stddev: pointer to the stddev value
 * todo: get result for 3 channels
 * *****************************************************************/
void meanstddev(float* input,int width, int height, int step1, float* mean, float* stddev)
{
    int number = width * height;
    float* data = (float *)input;
    float temp = 0;
    float sum = 0;
    int nn = width >>2;
    int nn_pow = nn<<2;
    int remain = width - nn_pow;
    int addr_remain = step1 - nn_pow;
    for(int i=0;i<height;i++)
    {
        nn = width >> 2;
        asm volatile(
             "mov w0,#0    \n"
             "dup v0.4s,w0 \n"
             "0:                     \n"
             "prfm pldl1keep,[%[src],#128]    \n"
             "ld1  {v1.4s},[%[src]]  \n"
             "fadd v0.4s,v0.4s,v1.4s \n"
             "add %[src],%[src],#16  \n"
             "subs %[i],%[i],#1    \n"
             "bne 0b               \n"
             "mov x0,v0.d[1] \n"
             "mov v1.d[0],x0 \n"
             "fadd v0.2s,v0.2s,v1.2s \n"
             "faddp v0.2s,v0.2s,v0.2s \n"
             "mov  %w[sum],v0.s[0] \n"
             :[i]"+r"(nn),
              [sum]"+r"(temp),
              [src]"+r"(data)
             :
             :"cc","memory","v0","v1","x0"
        );
        sum += temp;
        for(int j=0;j<remain;j++)
        {
            sum += input[nn_pow+j];
        }
        data += addr_remain;
    }
    float meanvalue = sum / number;
    temp = 0;
    sum =0;
    data =(float*)input;
    for(int i=0; i< height;i++)
    {
        nn = width >> 2;
        asm volatile(
             "mov w0,#0    \n"
             "dup v0.4s,w0 \n"
             "mov w0,%w[meanvalue] \n"
             "dup v2.4s,w0  \n"
             "0:                     \n"
             "prfm pldl1keep,[%[src],#128]    \n"
             "ld1  {v1.4s},[%[src]]  \n"
             "fsub v1.4s,v1.4s,v2.4s \n"
             "fmla v0.4s,v1.4s,v1.4s \n"
             "add %[src],%[src],#16  \n"
             "subs %[i],%[i],#1    \n"
             "bne 0b               \n"
             "mov x0,v0.d[1] \n"
             "mov v1.d[0],x0 \n"
             "fadd v0.2s,v0.2s,v1.2s \n"
             "faddp v0.2s,v0.2s,v0.2s \n"
             "mov  %w[sum],v0.s[0] \n"
             :[i]"+r"(nn),
              [sum]"+r"(temp),
              [src]"+r"(data)
             :[meanvalue]"r"(meanvalue)
             :"cc","memory","v0","v1","v2","x0"
        );
        sum += temp;
        for(int j=0;j<remain;j++)
        {
            sum += (input[nn_pow+j]-meanvalue)*(input[nn_pow+j]-meanvalue);
        }
        data += addr_remain;
    }
    float stdvalue = sqrt(sum/number);// the divisor should be n-1 (set to n to get the same result with opencv)
    *mean = meanvalue;
    *stddev = stdvalue;
}

void meanstddev_S16(short* input,int width, int height, int step1, float* mean, float* stddev)
{
    int number = width * height;
    short* data = (short *)input;
    int temp = 0;
    int sum = 0;
    int nn = width >>3;
    int nn_pow = nn<<3;
    int remain = width - nn_pow;
    int addr_remain = step1 - nn_pow;
    for(int i=0;i<height;i++)
    {
        nn = width >> 3;
        asm volatile(
             "mov w0,#0    \n"
             "dup v0.8h,w0 \n"
             "0:                     \n"
             "prfm pldl1keep,[%[src],#128]    \n"
             "ld1  {v1.8h},[%[src]]  \n"
             "add  v0.8h,v0.8h,v1.8h \n"
             "add %[src],%[src],#16  \n"
             "subs %[i],%[i],#1    \n"
             "bne 0b               \n"
             "saddlv s0,v0.8h \n"
             "fmov %w[sum],s0 \n"
             :[i]"+r"(nn),
              [sum]"+r"(temp),
              [src]"+r"(data)
             :
             :"cc","memory","v0","v1","x0","s0"
        );
        sum += temp;
        for(int j=0;j<remain;j++)
        {
            sum += input[nn_pow+j];
        }
        data += addr_remain;
    }
    float meanvalue = (float)sum / number;
    nn = width >>2;
    nn_pow = nn<<2;
    remain = width - nn_pow;
    addr_remain = step1 - nn_pow;
    float tempf = 0;
    float sumf =0;
    data =(short*)input;
    for(int i=0; i< height;i++)
    {
        nn = width >> 2;
        asm volatile(
             "mov w0,#0        \n"
             "dup v0.4s,w0     \n"
             "dup v3.4s,w0     \n"
             "dup v2.4s,%w[mean]     \n"
             "0:                     \n"
             "ld1  {v1.4h},[%[src]]  \n"
             "saddl  v1.4s,v1.4h,v3.4h \n" // !!! trick !!! 使用addl 加0的技巧把16位数转为32位  use movl get the same result
             "scvtf v1.4s,v1.4s \n"
             "fsub   v1.4s,v1.4s,v2.4s \n"
             "fmla   v0.4s,v1.4s,v1.4s \n"
             "add %[src],%[src],#8  \n"
             "subs %[i],%[i],#1    \n"
             "bne 0b               \n"
             "saddlv d0,v0.4s \n"
             "fmov %w[sum],s0 \n"
             :[i]"+r"(nn),
              [sum]"+r"(sumf),
              [src]"+r"(data),
              [mean]"+r"(meanvalue)
             :
             :"cc","memory","v0","v1","v2","v3","x0","d0"
        );
        sumf += tempf;
        for(int j=0;j<remain;j++)
        {
            tempf = (float)input[nn_pow +j]- meanvalue;
            sumf += tempf * tempf;
        }
        data += addr_remain;
    }
    float stdvalue = sqrt(sumf/number);// the divisor should be n-1 (set to n to get the same result with opencv)
    *mean = meanvalue;
    *stddev = stdvalue;
}


