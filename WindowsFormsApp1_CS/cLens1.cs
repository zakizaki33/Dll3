using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices; // *** これが必要

// namespace USING_Dll3_CS
namespace WindowsFormsApp1_CS
{
    public class cLens1
    {
        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern IntPtr Create_cLens1();
        
        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern void Delete_cLens1(IntPtr p_cLens1);
        
        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern void SetRadius(IntPtr p_cLens1, int i, double value);
        
        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern double GetRadius(IntPtr p_cLens1, int i);
        
        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern void SetDistance(IntPtr p_cLens1, int i, double value);
        
        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern double GetDistance(IntPtr p_cLens1, int i);
        
        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void SetGlassName(IntPtr p_cLens1, int i, string name);
        
        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        [return: MarshalAs(UnmanagedType.BStr)]
        public static extern string GetGlassName(IntPtr p_cLens1, int i);
        
        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern int GetK(IntPtr p_cLens1);
        
        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern double focallength(IntPtr p_cLens1);

        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern double backf(IntPtr p_cLens1);

        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern double frontf(IntPtr p_cLens1);

        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void SetStop(IntPtr p_cLens1, int i);

        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern int GetStop(IntPtr p_cLens1);

        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern void SetEAy(IntPtr p_cLens1, int i, double value);

        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern double GetEAy(IntPtr p_cLens1, int i);

        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern void EPCalculation(IntPtr p_cLens1);

        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern void Set_s(IntPtr p_cLens1, double value);

        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern double Get_s(IntPtr p_cLens1);

        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern void Set_t(IntPtr p_cLens1, double value);

        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern double Get_t(IntPtr p_cLens1);

        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern void Set_EPD(IntPtr p_cLens1, double value);

        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern double Get_EPD(IntPtr p_cLens1);

        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern void SetColorN(IntPtr p_cLens1, int num);

        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern int GetColorN(IntPtr p_cLens1);

        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern void SetColorWeight(IntPtr p_cLens1, int i, double value);

        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern double GetColorWeight(IntPtr p_cLens1, int i);

        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void SetColor(IntPtr p_cLens1, int i, string name);

        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        [return: MarshalAs(UnmanagedType.BStr)]
        public static extern string GetColor(IntPtr p_cLens1, int i);

        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern void MakeSAGraph(IntPtr p_cLens1);

        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern void SaveAsBmp(IntPtr p_cLens1);

        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern void MakeLensView(IntPtr p_cLens1);

        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern void SaveAsBmpLensView(IntPtr p_cLens1);

        // DLL間でやり取りをするためのポインタを定義
        public IntPtr _cLens1Pointer;

        public cLens1()
        {
            _cLens1Pointer = Create_cLens1();
        }

        ~cLens1()
        {
            Delete_cLens1(_cLens1Pointer);
        }

        // 名前をC++と同じにして良い(成立する)
        public void SetRadius(int i, double value)
        {
            SetRadius(_cLens1Pointer, i, value);
        }
        public double GetRadius(int i)
        {
            return GetRadius(_cLens1Pointer, i);
        }
        public void SetDistance(int i, double value)
        {
            SetDistance(_cLens1Pointer, i, value);
        }
        public double GetDistance(int i)
        {
            return GetDistance(_cLens1Pointer, i);
        }
        public void SetGlassName(int i, string name)
        {
            SetGlassName(_cLens1Pointer, i, name);
        }
        public string GetGlassName(int i)
        {
            return GetGlassName(_cLens1Pointer, i);
        }
        public int GetK()
        {
            return GetK(_cLens1Pointer);
        }
        public double focallength()
        {
            return focallength(_cLens1Pointer);
        }
        public double backf()
        {
            return backf(_cLens1Pointer);
        }
        public double frontf()
        {
            return frontf(_cLens1Pointer);
        }
        public void SetStop(int i)
        {
            SetStop(_cLens1Pointer, i);
        }
        public int GetStop()
        {
            return GetStop(_cLens1Pointer);
        }
        public void SetEAy(int i, double value)
        {
            SetEAy(_cLens1Pointer, i, value);
        }
        public double GetEAy(int i)
        {
            return GetEAy(_cLens1Pointer, i);
        }
        public void EPCalculation()
        {
            EPCalculation(_cLens1Pointer);
        }
        public void Set_s(double value)
        {
            Set_s(_cLens1Pointer, value);
        }
        public double Get_s()
        {
            return Get_s(_cLens1Pointer);
        }
        public void Set_t(double value)
        {
            Set_t(_cLens1Pointer, value);
        }
        public double Get_t()
        {
            return Get_t(_cLens1Pointer);
        }
        public void Set_EPD(double value)
        {
            Set_EPD(_cLens1Pointer, value);
        }
        public double Get_EPD()
        {
            return Get_EPD(_cLens1Pointer);
        }
        public void SetColorN(int i)
        {
            SetColorN(_cLens1Pointer, i);
        }
        public int GetColorN()
        {
            return GetColorN(_cLens1Pointer);
        }
        public void SetColorWeight(int i, double value)
        {
            SetColorWeight(_cLens1Pointer, i, value);
        }
        public double GetColorWeight(int i)
        {
            return GetColorWeight(_cLens1Pointer, i);
        }
        public void SetColor(int i, string name)
        {
            SetColor(_cLens1Pointer, i, name);
        }
        public string GetColor(int i)
        {
            return GetColor(_cLens1Pointer, i);
        }
        public void MakeSAGraph()
        {
            MakeSAGraph(_cLens1Pointer);
        }
        public void SaveAsBmp()
        {
            SaveAsBmp(_cLens1Pointer);
        }
        public void MakeLensView()
        {
            MakeLensView(_cLens1Pointer);
        }
        public void SaveAsBmpLensView()
        {
            SaveAsBmpLensView(_cLens1Pointer);
        }
    }
}
