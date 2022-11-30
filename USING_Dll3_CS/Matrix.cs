using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices; // *** これが必要

namespace USING_Dll3_CS
{
    public class Matrix2
    {
    
        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern IntPtr CreateMatrix();
        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern IntPtr CreateMatrixMN(int m, int n);
        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern void SetMatrixIJ(IntPtr pMatrix, int m, int n, double value);
        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern double GetMatrixIJ(IntPtr pMatrix, int m, int n);
        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern void DeleteMatrix(IntPtr pMatrix);

        private readonly IntPtr _MatrixPointer;

        public Matrix2()
        {
           _MatrixPointer = CreateMatrix();
        }

        public Matrix2(int m, int n)
        {
            _MatrixPointer = CreateMatrixMN(m, n);
        }

        ~Matrix2()
        {
            DeleteMatrix(_MatrixPointer);
        }

        public void SetMatrix(int m, int n, double value)
        {
            SetMatrixIJ(_MatrixPointer, m, n, value);
        }
        public double GetMatrix(int m, int n)
        {
            return GetMatrixIJ(_MatrixPointer, m, n);
        }
    }
}
