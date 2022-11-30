using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices; // *** これが必要

namespace USING_Dll3_CS
{
    class Program
    {
        [StructLayout(LayoutKind.Sequential)]
        public struct Complex 
        {
            public double x;
            public double y;
            // public IntPtr data;
        }

        [StructLayout(LayoutKind.Sequential)]
        public struct Matrix
        {
            public int n;
            public int m;
            public IntPtr data;
        }

        // 追加関数
        // [DllImport("Dll3.dll")]
        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void test_01();
        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern double ARG(ref Complex pComp);
        [DllImport("Dll_CPPtoCS.dll")]
        public static extern void func();
        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void Conjugate(ref Complex cComp1, ref Complex cComp2);
        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void matrix_init(ref Matrix p1, int m, int n);
        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void matrix_inv(IntPtr p1, IntPtr p2);



        static void Main(string[] args)
        {

            Console.WriteLine("最初の一歩");
            func(); // Dll_CPPtoCS.dll の中の関数
            test_01();　// Dll3.dll の中の関数

            // complexを計算させる
            var comp1 = new Complex(){ x = Math.Sqrt(3), y = 1 };
            Console.WriteLine($"{comp1.x} {comp1.y}");

            double arg = ARG(ref comp1);
            Console.WriteLine(arg*180/Math.PI);

            // 複素共役を取得する（インプットとアウトプットの二つの複素数で戻り値無し）
            // 空の複素数を用意
            Complex comp_conj = new Complex() {x=0, y=0};
            Conjugate(ref comp1, ref comp_conj);
            Console.WriteLine("複素共役をゲットする");
            Console.WriteLine("comp1実部=" + $"{comp1.x}" + " ,虚部=" + $"{comp1.y}\n");
            Console.WriteLine("comp1共役実部=" + $"{comp_conj.x}" + " ,comp1共役虚部=" + $"{comp_conj.y}\n");


            // Matrix の関数を利用する
            // var c4 = new Matrix();
            // matrix_init(ref c4, 2, 2);

            Matrix2 matrix1 = new Matrix2();

            Matrix2 matrix2 = new Matrix2(2,2);

            // 行列に値をセット
            matrix2.SetMatrix(1,1,2);
            matrix2.SetMatrix(1, 2, 1);
            matrix2.SetMatrix(2, 1, 7);
            matrix2.SetMatrix(2, 2, 4);

            Console.WriteLine(matrix2.GetMatrix(1,1));

            // コピー用のMatrixを作る
            Matrix2 matrix2_copy = new Matrix2(2, 2);

            matrix_inv(matrix2._MatrixPointer, matrix2_copy._MatrixPointer);

            Console.WriteLine(matrix2_copy.GetMatrix(1, 1));
            Console.WriteLine(matrix2_copy.GetMatrix(1, 2));
            Console.WriteLine(matrix2_copy.GetMatrix(2, 1));
            Console.WriteLine(matrix2_copy.GetMatrix(2, 2));

        }
    }
}
