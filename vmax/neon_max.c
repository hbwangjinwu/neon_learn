#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//v8支持maxv 指令直接求出向量寄存器中的最大值，而在v7 架构中需要借助通用寄存器复制 然后使用maxp 指令求出最大值
int main()
{
    //neon汇编不支持使用数组名作为指针传递，汇编中会改变num的值
    float num[16] = {1,3,4,5,6,2,7,7,4,10,12,13,15,16,18};
    //float* src = (float*)malloc(sizeof(float)* 16);
    //float* src_old = src;
    //memcpy(src,num,sizeof(float)*16);
    int nn = 4;
    float max = 0;
    asm volatile(
            "mov w0,#0      \n"
            "dup v0.4s,w0   \n"//set v0 to 0
            "0:             \n"
            //"ld1 {v1.4s},[%[src]] \n"
            "ld1 {v1.4s},[%[num]], x0 \n"
            "add x0 , x0 ,#16 \n"
            "fmax v0.4s,v0.4s,v1.4s \n" //fmax or smax max
            //"add  %w[src],%w[src],#16 \n"
            "subs %w[i],%w[i],#1 \n" //subs 影响状态寄存器
            "bne  0b             \n" //需要使用0b 进行跳转
            "fmaxv s0,v0.4s      \n" //use fmaxv/smaxv get max value cross the vector register v8 only
            "fmov  %w[max],s0    \n" //usr fmov from s0 to reguler register v8 only
            :[max]"+r"(max),
             [i]"+r"(nn)
             //[src]"+r"(src) //+read write | & use only for output |  = write only | r tell gcc to  specify a common register
            :[num]"r"(num) //如果不使用malloc 分配内存的话可以把num传递到这里， ld1的时候记得加上offset
            :"cc","memory","v0","v1","w0","s0"
    );
    printf("max value is %f \n",max);
    //free(src_old);//这里如果free(src) 会出现double free 的问题
}


