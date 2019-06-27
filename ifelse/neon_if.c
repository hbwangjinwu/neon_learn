#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main()
{
    float num[16] = {1,-3,4,5,-6,2,7,-7,4,-10,12,-13,15,-16,18,20};
    int nn = 4;
    float *addr = num;
    asm volatile(
            "mov w0,#0      \n"
            "dup v0.4s,w0   \n"//set v0 to 0
            "0:             \n"
            "ld1 {v1.4s},[%[num]]    \n"
            "cmgt v2.4s,v1.4s,v0.4s  \n"
            "bsl  v2.16b,v1.16b,v0.16b  \n"//bsl must byte
            "st1 {v2.4s},[%[num]],#16 \n"
            "subs %w[i],%w[i],#1 \n" //subs 影响状态寄存器
            "bne  0b             \n" //需要使用0b 进行跳转
            :[i]"+r"(nn),
             [num]"+r"(addr)
            :
            :"cc","memory","v0","v1","v2","w0"
    );
    for(int i=0;i<16;i++)
    {
       printf("%f \n",num[i]);
    }
}
