/* Many people have requested a simple example on how to create a C
 * MEX-file.  In response to this request, the following C MEX-file,
 * named mexample, is provided as an introduction to cmex
 * programming. mexample is a commented program which describes how to
 * use the following MEX-functions:
 *
 *许多人要求提供一个简单的例子来说明如何创建c mex文件。为了响应此请求，以下名为
 * mexample的c mex文件作为c mex编程的介绍提供。mexample是一个注释程序，描述如
 * 何使用以下mex函数：
 *
 * mexErrMsgTxt
 * mxCreateDoubleMatrix
 * mxGetM
 * mxGetN
 * mxGetPr
 * mxIsComplex
 * mxIsSparse
 * mxIsChar
 *
 * In MATLAB, mexample accepts two inputs and returns one output. The
 * inputs are a 2x2 array denoted as ARRAY_IN and a 2x1 vector denoted as
 * VECTOR_IN.  The function calculates the determinant of ARRAY_IN,
 * multiplies each element of VECTOR_IN by the determinant, and returns
 * this as the output, denoted by VECTOR_OUT.  All inputs and outputs to
 * this function are assumed to be real (not complex). 
 *在Matlab中，Mexample接受两个输入并返回一个输出。输入是表示为array_in的2x2数
 * 组和表示为vector_in的2x1向量。该函数计算数组的行列式，将向量的每个元素乘以
 * 行列式，并将其作为输出返回，用向量表示。此函数的所有输入和输出都假定为真实的
 * （不复杂）。
 */

/*   First, include some basic header files.  The header file
 * "mex.h" is required for a MEX-file.  Add any other header
 * files that your function may need here.
 *首先，包括一些基本的头文件。mex文件需要头文件“mex.h”。在此处添加您的函数可
 * 能需要的任何其他头文件。
 */

#include "mex.h"
#include <limits.h>
#include <float.h>
#include <math.h>       /* log */
/*   A C MEX-file generally consists of two sections.  The first
 * section is a function or set of functions which performs
 * the actual mathematical calculation that the MEX-function
 * is to carry out.  In this example, the function is called
 * workFcn().  The second section is a gateway between MATLAB
 * and the first section, and consists of a function called
 * mexFunction.  The gateway is responsible for several tasks,
 * including:
 *c mex文件通常由两部分组成。第一部分是一个函数或一组函数，用于执行MEX函数要
 * 执行的实际数学计算。在本例中，函数名为workfcn（）。第二部分是Matlab和第一
 * 部分之间的网关，由一个名为mexFunction的函数组成。网关负责几个任务，包括：
 *
 * I)  error checking,
 * II)  allocating memory for return arguments,
 * III)  converting data from MATLAB into a format that
 * the workFcn function can use, and vice versa.
 *i）错误检查，
 * ii）为返回参数分配内存，
 * iii）将来自Matlab的数据转换为WorkFCN函数可以使用的格式，反之亦然。

 *
 * The first function to be written in this example, then, is
 * workFcn:
 *本例中要编写的第一个函数是workfcn：
 *
 * Since C and MATLAB handle two-dimensional arrays
 * differently, we will explicitly declare the dimension of
 * the variable theArray.  The variables, theVector and
 * theResult, are both one-dimensional arrays, and therefore
 * do not need such rigid typing. 
 *由于C和Matlab处理二维数组的方式不同，我们将显式地声明数组变量的维数。
 * 变量vector和result都是一维数组，因此不需要如此严格的类型。
 */


void viterbi(
        int N,
        int T,
        double a_matrix[4][4],
        int max_duration_D,
        double *delta,
        double *observation_probs,
        double duration_probs [4][150],
        double *psi,
        double *psi_duration_out,
        double duration_sum_in[4]
        )
        
{
    
    int i;
    int i2;
    int i3;
    int j;
    int t;
    int d;
    
    
    
    for (t = 1; t<T+ max_duration_D-1;t++){
        
        
        /*For each state 对于每个状态*/
        for (j = 0; j<4;j++){
            
            double emission_probs = 0;
            
            /*        max_duration_D*/
            for (d = 1; d<=max_duration_D; d++){
                
                int start; int max_index = 0;
                int end_t = 0;
                float probs = 0;
                float duration_sum = 0;
                float delta_temp = 0;
                float max_delta = -1*DBL_MAX;
                
                
                /*  Get the maximum value for delta at this t, and record the state where it was found:
                 * This is the first half of the expression of equation 33a from Rabiner:
                 在t处获得delta的最大值，并记录它被找到的状态：这是Rabiner公式表达式的前半部分。*/
                
                /*
                 * %The start of the analysis window, which is the current time
                 * %step, minus d (the time horizon we are currently looking back),
                 * %plus 1. The analysis window can be seen to be starting one
                 * %step back each time the variable d is increased.
                 * % This is clamped to 1 if extending past the start of the
                 * % record, and T-1 is extending past the end of the record:
                 *分析窗口的开始，即当前时间点，减去d（我们当前回顾的时间范围），再加1。
                 * 可以看到分析窗口在每次增加变量d时都会向后退一步。如果延伸超过记录的开始，
                 * 则将其钳制为1，而t-1延伸超过记录的结束：
                 */
                
                start = t - d;
                
                if(start < 0){
                    start = 0;
                }
                
                if(start > T-2){
                    start = T-2;
                }
                
                /*
                 * %The end of the analysis window, which is the current time
                 * %step, unless the time has gone past T, the end of the record, in
                 * %which case it is truncated to T. This allows the analysis
                 * %window to extend past the end of the record, so that the
                 * %timing durations of the states do not have to "end" at the end
                 * %of the record.
                 *分析窗口的结尾，即当前时间步，除非时间已超过t，否则为记录的结尾，在这种情况下，
                 * 它将被截断为t。这允许分析窗口扩展到记录的结尾，这样状态的计时持续时间就不必
                 * “结束”在记录的结尾。
                 */
                
                end_t = t;
                if(end_t>T-1){
                    end_t = T-1;
                }
                
                
                for(i = 0; i<N; i++)
                {
                    double temp = delta[(start) +(i*(T+ max_duration_D-1))] + log(a_matrix[i][j]);
                    if(temp > max_delta){
                        max_delta = temp;
                        max_index = i;
                    }
                }
                
                
                /*//Find the normaliser for the observations at the start of the
                 * //analysis window. The probability of seeing all the
                 * //observations in the analysis window in state j is updated each
                 * //time d is incrememented two lines below, so we only need to
                 * //find the observation probabilities for one time step, each
                 * //time d is updated:
                 *在分析窗口的开始处找到观测值的规范器。在分析窗口中看到状态j中的所有观测值
                 *的概率在每次d增加以下两行时更新，因此我们只需要找到一个时间步的观测概率，每次d都更新：
                 */
                
                
                probs = 0;
                for(i2 = start; i2<=end_t; i2++){
                    
                    // Ensure that the probabilities aren't zero leading to -inf probabilities after log:
                    //确保概率不为零，从而在日志之后导致-inf概率
                    if(observation_probs[i2 +j*T] == 0){
                        observation_probs[i2 +j*T] = FLT_MIN;
                    }
                    
                    probs = probs + log(observation_probs[i2 +j*T]);
                }
                
                if(probs ==0){
                    probs = FLT_MIN;
                }

                emission_probs = (probs);
                
                /*Find the total probability of transitioning from the last
                 * //state to this one, with the observations and being in the same
                 * //state for the analysis window. This is the duration-dependant
                 * //variation of equation 33a from Rabiner:
                 *找到从上一个状态转换到这个状态的总概率，观察结果与分析窗口处于相同的状态。
                 * 这是rabiner方程33a的持续时间相关变化
                 */
                delta_temp = max_delta + (emission_probs)+ (log((duration_probs[j][d-1]/duration_sum_in[j])));
                
                
                
                // Uncomment the below for debuggin:为debuggin取消下面的注释:
//                     mexPrintf("\n t:%d", t);
//                     mexPrintf("\n j:%d", j);
//                     mexPrintf("\n d:%d", d);
//                     mexPrintf("\n max_delta:%f", max_delta);
//                     mexPrintf("\n max_index:%i \n", max_index);
//                     mexPrintf ("emission_probs: %f \n",emission_probs);
//                     mexPrintf ("log((duration_probs[j][d-1]/duration_sum)): %f \n",log((duration_probs[j][d-1]/duration_sum_in[j])));
//                     mexPrintf ("delta_temp: %f \n",delta_temp);
//                     mexPrintf ("delta[t+j*(T+ max_duration_D-1)]: %f \n",delta[t+j*(T+ max_duration_D-1)]);
//                     mexPrintf ("duration_probs[j][d]: %f \n",duration_probs[j][d]);
//                     mexPrintf ("duration_sum_in[j]: %f \n",duration_sum_in[j]);
                  
                /*
                 * Unlike equation 33a from Rabiner, the maximum delta could come
                 * from multiple d values, or from multiple size of the analysis
                 * window. Therefore, only keep the maximum delta value over the
                 * entire analysis window:
                 * If this probability is greater than the last greatest,
                 * update the delta matrix and the time duration variable:
                 *与Rabiner方程33a不同，最大增量可以来自多个D值，或者来自分析窗口的多个大小。
                 *因此，只将最大增量值保留在整个分析窗口上：如果此概率大于最后一个最大值，
                 *则更新delta矩阵和持续时间变量：
                 */
                
                if(delta_temp>delta[t+j*(T+ max_duration_D-1)]){
                    
                    delta[t+j*(T+ max_duration_D-1)] = delta_temp;
                    psi[t+j*(T+ max_duration_D-1)] = max_index+1;
                    
                    psi_duration_out[t + j*(T+ max_duration_D-1)] = d;
                    
                }
            }
        }
    }
    
}

/*   Now, define the gateway function, i.e., mexFunction.Below
 * is the standard, predeclared header to mexFunction.  nlhs
 * and nrhs are the number of left-hand and right-hand side
 * arguments that mexample was called with from within MATLAB.
 * In this example, nlhs equals 1 and nrhs should equal 2.  If
 * not, then the user has called mexample the wrong way and
 * should be informed of this.  plhs and prhs are arrays which
 * contain the pointers to the MATLAB arrays, which are
 * stored in a C struct called an Array.  prhs is an array of
 * length rhs,and its pointers point to valid input data.
 * plhs is an array of length nlhs, and its pointers point to
 * invalid data (i.e., garbage).  It is the job of mexFunction
 * to fill plhs with valid data.
 *现在，定义网关功能，即mexfunction。下面是mexfunction的标准预分离头。
 * nlhs和nrhs是从Matlab中调用Mexample时使用的左侧和右侧参数的数目。
 * 在本例中，nlhs等于1，nrhs应等于2。如果没有，则用户以错误的方式调用了mexample，
 * 应将此情况告知用户。plh和prh是包含指向matlab数组的指针的数组，这些指针存储在一个称为数组的c结构中。
 * prhs是一个长度为rhs的数组，它的指针指向有效的输入数据。plhs是长度nlhs的数组，它的指针指向无效数据（即垃圾）。
 * 用有效数据填充PLH是MexFunction的工作。
 *
 * First, define the following values.  This makes it much
 * easier to change the order of inputs to mexample, should we
 * want to change the function later.  In addition, it makes
 * the code easier to read.
 *首先，定义以下值。如果我们想稍后更改函数，这使得将输入顺序更改为mexample变得容易得多。
 * 此外，它使代码更易于阅读。
 */

#define N prhs[0]
#define T prhs[1]
#define a_matrix prhs[2]
#define max_duration_D prhs[3]
#define delta  prhs[4]
#define observation_probs prhs[5]
#define duration_probs prhs[6]
#define psi prhs[7]
#define duration_sum prhs[8]


#define delta_out plhs[0]
#define psi_out plhs[1]
#define psi_duration plhs[2]


void mexFunction(
        int     nlhs,
        mxArray  *plhs[],
        int     nrhs,
        const mxArray  *prhs[]
        )
{
    double a_matrix_in[4][4];/* 2 dimensional C array to pass to workFcn() */
    double *delta_in_matrix;/* 2 dimensional C array to pass to workFcn() */
    double *observation_probs_matrix;/* 2 dimensional C array to pass to workFcn() */
    double *psi_matrix;/* 2 dimensional C array to pass to workFcn() */
    double duration_sum_in[4];/* 2 dimensional C array to pass to workFcn() */
    
    double duration_probs_matrix[4][150];/* 2 dimensional C array to pass to workFcn() */
    
    int actual_T;
    int fake_T_extended;
    int actual_N;
    int max_duration_D_val;
    
    int    row,col;        /* loop indices */
    int    m,n;            /* temporary array size holders */
    
    /*   Step 1: Error Checking Step 1a: is nlhs 1?  If not,
     * generate an error message and exit mexample (mexErrMsgTxt
     * does this for us!)
     *步骤1：错误检查步骤1a：是nlhs 1吗？如果不是，生成一个错误消息并退出示例（MexErrsGtxt为我们这样做！）
     */
    if (nlhs!=3)
        mexErrMsgTxt("mexample requires 3 output argument.");
    
    /*   Step 1b: is nrhs 2? 步骤1b是nrhs2吗？*/
    if (nrhs!=9)
        mexErrMsgTxt("mexample requires 9 input arguments");
    
    
    actual_T = mxGetM(observation_probs);
    actual_N = mxGetN(observation_probs);
    
    max_duration_D_val = mxGetScalar(max_duration_D);
    
    
    /*   Step 2:  Allocate memory for return argument(s) 步骤2：为返回参数分配内存*/
    delta_out = mxCreateDoubleMatrix((actual_T+max_duration_D_val-1), actual_N, mxREAL);
    psi_out = mxCreateDoubleMatrix((actual_T+max_duration_D_val-1), actual_N, mxREAL);
    psi_duration = mxCreateDoubleMatrix((actual_T+max_duration_D_val-1), actual_N, mxREAL);
    
    /*   Step 3:  Convert ARRAY_IN to a 2x2 C array
     * MATLAB stores a two-dimensional matrix in memory as a one-
     * dimensional array.  If the matrix is size MxN, then the
     * first M elements of the one-dimensional array correspond to
     * the first column of the matrix, and the next M elements
     * correspond to the second column, etc. The following loop
     * converts from MATLAB format to C format:
     *第三步：将数组转换成2x2c数组，Matlab将二维矩阵作为一维数组存储在内存中。
     * 如果矩阵是大小mxn，则一维数组的前m个元素对应于矩阵的第一列，下m个元素对应于第二列等。
     * 以下循环将从Matlab格式转换为C格式：
     */
    
    for (col=0; col < mxGetN(a_matrix); col++){
        for (row=0; row < mxGetM(a_matrix); row++){
            a_matrix_in[row][col] =(mxGetPr(a_matrix))[row+col*mxGetM(a_matrix)];
        }
    }
    
    for (col=0; col < mxGetM(duration_sum); col++){
        duration_sum_in[col] =(mxGetPr(duration_sum))[col];
    }
    
    
    
    
    delta_in_matrix = mxGetPr(delta);
    observation_probs_matrix = mxGetPr(observation_probs);
    psi_matrix = mxGetPr(psi);
    
    /*     for (col=0; col < mxGetN(delta); col++){
     * //         for (row=0; row < mxGetM(delta); row++){
     * //
     * //
     * //             observation_probs_matrix[row][col] =(mxGetPr(observation_probs))[row+col*mxGetM(observation_probs)];
     * //             psi_matrix[row][col] =(mxGetPr(psi))[row+col*mxGetM(psi)];
     * //         }
     * //     }*/
    
    
    for (col=0; col < mxGetN(duration_probs); col++){
        for (row=0; row < mxGetM(duration_probs); row++){
            duration_probs_matrix[row][col] =(mxGetPr(duration_probs))[row+col*mxGetM(duration_probs)];
        }
    }
    
    
    
    
    /*   mxGetPr returns a pointer to the real part of the array
     * ARRAY_IN.  In the line above, it is treated as the one-
     * dimensional array mentioned in the previous comment. 
     *mxgetpr返回指向数组数组实际部分的指针。在上面的行中，它被视为前面注释中提到的一维数组。
     */
    
    /*   Step 4:  Call workFcn function 步骤4：调用workfcn函数*/
    viterbi(actual_N,actual_T,a_matrix_in,max_duration_D_val,delta_in_matrix,observation_probs_matrix,duration_probs_matrix,psi_matrix,mxGetPr(psi_duration),duration_sum_in);
    memcpy ( mxGetPr(delta_out), delta_in_matrix, actual_N*(actual_T+max_duration_D_val-1)*8);
    memcpy ( mxGetPr(psi_out), psi_matrix, actual_N*(actual_T+max_duration_D_val-1)*8);
    
}