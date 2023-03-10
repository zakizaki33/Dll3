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

        public IntPtr _cLens1Pointer;

        public cLens1()
        {
            _cLens1Pointer = Create_cLens1();
        }

        ~cLens1()
        {
            Delete_cLens1(_cLens1Pointer);
        }

        // 名前をC++と同じにして良いのだろうか。。。。
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
    }
}
