﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices; // *** これが必要

namespace USING_Dll3_CS
{
    class Program
    {
        
        

        // 追加関数
        // [DllImport("Dll3.dll")]
        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void test_01();
        [DllImport("Dll_CPPtoCS.dll")]
        public static extern void func();


        static void Main(string[] args)
        {

            Console.WriteLine("最初の一歩");
            func(); // Dll_CPPtoCS.dll の中の関数
            test_01();　// Dll3.dll の中の関数

        }
    }
}