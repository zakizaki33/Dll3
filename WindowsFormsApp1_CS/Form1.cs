using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Xml.Serialization;
using System.Windows;
using System.Collections.ObjectModel;
using Microsoft.VisualBasic.FileIO;
using System.Runtime.InteropServices; // *** これが必要

namespace WindowsFormsApp1_CS
{
    public partial class Form1 : Form
    {
        // 変数の設定
        int x_Rect = 111;
        int y_Rect =35;
        
        // https://www.wareko.jp/blog/output-text-string-to-console-window-with-windows-form-application-in-c-sharp
        [System.Runtime.InteropServices.DllImport("kernel32.dll")] // この行を追加
        private static extern bool AllocConsole();                 // この行を追加
        // [System.Runtime.InteropServices.DllImport("MyDllOptics.dll")] // この行を追加
        // private static extern bool AllocConsole();                 // この行を追加

        // 追加関数 2023-03-02
        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void test_01();
        [DllImport("Dll3.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern double Return123();


        // Lensクラスインスタンス生成
        // var lens = new ClassLibrary1.Lens();
        readonly ClassLibrary1.Lens lens = new ClassLibrary1.Lens();
        // readonly MydllOptics.Lens lens2 = new MydllOptics.Lens();

        // cLens1の生成のテスト　
        cLens1 plens1 = new cLens1();
        // Userクラスをまとめて処理する準備
        ObservableCollection<User> data = new ObservableCollection<User>();

        public Form1()
        {
            InitializeComponent();
            AllocConsole();

            // コンソールに出力されることを確認する
            Console.WriteLine("コンソールへの出力を確認！！！　2023-03-02");

            // Dll3の関数の呼び出し
            System.Diagnostics.Debug.WriteLine("テストメッセージ");
            Console.WriteLine("test_01");
            test_01();　// Dll3.dll の中の関数　（なぜか表示されない。。。。）
            Console.WriteLine(Return123());
            
            /*
            // 値のセット
            plens1.SetRadius(1, 100);
            plens1.SetRadius(2, -100);
            plens1.SetDistance(1, 10);
            // 値の確認
            Console.WriteLine("第1面の値を確認 ⇒ " + $"{plens1.GetRadius(1)}\n");
            Console.WriteLine("第2面の値を確認 ⇒ " + $"{plens1.GetRadius(2)}\n");
            Console.WriteLine("第1面の厚さを確認 ⇒ " + $"{plens1.GetDistance(1)}\n");
            Console.WriteLine("面の総数を確認 ⇒ " + $"{plens1.GetK()}\n");


            string name1 = "518640";
            // string name1 = "500640";

            plens1.SetGlassName(1, name1);
            Console.WriteLine("第1面のガラス名を確認 ⇒ " + plens1.GetGlassName(1)+"だによ");
            
            Console.WriteLine("焦点距離を確認 ⇒ " + $"{plens1.focallength()}\n");
            */

            // RDNをセットしていく
            // RD は1から開始
            // Nはゼロから開始
            lens.k = 2; // 面の数をセットする

            lens.set_r(1, 100);
            lens.set_r(2, -100);
            lens.set_d(1, 10);
            lens.set_N(0, 1);
            lens.set_N(1, 1.5);
            lens.set_N(2, 1);

            // ドロップを許可
            pictureBox1.AllowDrop = true;

            // PictureBox2にお絵かきをする
            Line_Something();

            // dataの練習(Userクラスをまとめての取り扱い)
            data.Add(new User("太郎", "北海道"));
            data.Add(new User("次\"郎", "岩手"));
            data.Add(new User("三郎", "宮城"));
            // dataGridView1.DataSource = data;  //ここの工程で、datagridviewに書き込まれる
        }

        private void Form1_Load(object sender, EventArgs e)
        {
            // https://www.sejuku.net/blog/58436
            // カラム数を指定
            dataGridView1.ColumnCount = 4;

            // カラム名を指定
            dataGridView1.Columns[0].HeaderText = "Radius";
            dataGridView1.Columns[1].HeaderText = "Distance";
            dataGridView1.Columns[2].HeaderText = "Nd";
            dataGridView1.Columns[3].HeaderText = "Aperture";

            // データを追加(まずは初期値)
            dataGridView1.Rows.Add(100, 0, "1", 2.0);
            dataGridView1.Rows.Add(-100, 10, "518640", 2.0);
            dataGridView1.Rows.Add(0, 0, "1", 2.0);

            // textboxの初期化
            textBox1.Text = "練習 (^^♪";

            // Form初期化時によばれるもの、Formをダブルクリックすると生成できる
            // 矩形描画
            // pictureBox2.Paint += pictureBox2_Paint; // 主導でイベントを設定する場合、Formのデザインでイベント側で追加してもよい


        }

        private void button1_Click(object sender, EventArgs e)
        {
            // cLens1のdllを使って焦点距離が計算できるようにする
            // レンズデータを読み込む
            /*
            lens.set_r(1, Convert.ToDouble(dataGridView1.Rows[0].Cells[0].Value));
            lens.set_r(2, Convert.ToDouble(dataGridView1.Rows[1].Cells[0].Value));
            lens.set_d(1, Convert.ToDouble(dataGridView1.Rows[1].Cells[1].Value));
            lens.set_N(0, Convert.ToDouble(dataGridView1.Rows[0].Cells[2].Value));
            lens.set_N(1, Convert.ToDouble(dataGridView1.Rows[1].Cells[2].Value));
            lens.set_N(2, Convert.ToDouble(dataGridView1.Rows[2].Cells[2].Value));
            
            System.Console.WriteLine(lens.FocalLength());   // lensをどこでインスタンス化させてホールドし続ける？？
            label1.Text = "Focla Length =" + lens.FocalLength();
             */
            plens1.SetRadius(1, Convert.ToDouble(dataGridView1.Rows[0].Cells[0].Value));
            plens1.SetRadius(2, Convert.ToDouble(dataGridView1.Rows[1].Cells[0].Value));
            plens1.SetDistance(1, Convert.ToDouble(dataGridView1.Rows[1].Cells[1].Value));
            plens1.SetGlassName(0, Convert.ToString(dataGridView1.Rows[0].Cells[2].Value));
            plens1.SetGlassName(1, Convert.ToString(dataGridView1.Rows[1].Cells[2].Value));
            plens1.SetGlassName(2, Convert.ToString(dataGridView1.Rows[2].Cells[2].Value));
            plens1.SetEAy(1, Convert.ToDouble(dataGridView1.Rows[0].Cells[3].Value));
            plens1.SetEAy(2, Convert.ToDouble(dataGridView1.Rows[1].Cells[3].Value));

            // EPD等がただしく計算出来ているか確認
            plens1.SetStop(2);
            System.Console.WriteLine("絞り面は⇒" + $"{plens1.GetStop()}\n");
            plens1.EPCalculation();
            System.Console.WriteLine("s(物体面距離)は⇒" + $"{plens1.Get_s()}\n");
            System.Console.WriteLine("t(入射瞳位置)は⇒" + $"{plens1.Get_t()}\n");
            System.Console.WriteLine("EPDは⇒" + $"{plens1.Get_EPD()}\n");
            // 波長関連動作確認
            plens1.SetColor(1, "d"); // 一旦d線に変更
            System.Console.WriteLine("波長の数は⇒" + $"{plens1.GetColorN()}\n");
            System.Console.WriteLine("波長名前は⇒" + $"{plens1.GetColor(1)}\n");
            System.Console.WriteLine("波長の重みは⇒" + $"{plens1.GetColorWeight(1)}\n");
            plens1.SetColor(1, "e"); // デフォルトのe線に戻す
            // 焦点距離確認
            System.Console.WriteLine("焦点距離⇒" + $"{plens1.focallength()}\n");  
            label1.Text = "Focal Length =" + plens1.focallength() + "\n"+
                          "Back Focal Length =" + plens1.backf() + "\n";

            // 画像出力確認
            plens1.MakeSAGraph();
            plens1.SaveAsBmp();
            // pictureboxへ表示させる
            pictureBox1.SizeMode = PictureBoxSizeMode.Zoom;
            pictureBox1.ImageLocation = @"C:\Users\13273_Yamazaki\source\repos\Dll3\WindowsFormsApp1_CS\bin\x64\Debug\SAGraph_test_dll3.bmp";
        }

        private void dataGridView1_CellContentClick(object sender, DataGridViewCellEventArgs e)
        {

        }

        private void label1_Click(object sender, EventArgs e)
        {

        }

        private void textBox1_TextChanged(object sender, EventArgs e)
        {

        }

        private void pictureBox1_DragDrop(object sender, DragEventArgs e)
        {
            // https://chigusa-web.com/blog/cs-picturebox-dd/
            // ファイルパスの取得
            string[] fileName = (string[])e.Data.GetData(DataFormats.FileDrop);
            // 複数ファイルの場合最初の一個を表示
            pictureBox1.ImageLocation = fileName[0];
        }

        private void pictureBox1_DragEnter(object sender, DragEventArgs e)
        {
            if (e.Data.GetDataPresent(DataFormats.FileDrop))
            {
                e.Effect = DragDropEffects.All;
            }
            else
            {
                e.Effect = DragDropEffects.None;
            }
        }

        private void dataGridView1_DragDrop(object sender, DragEventArgs e)
        {

        }

        private void dataGridView1_DragEnter(object sender, DragEventArgs e)
        {

        }

        private void button2_Click(object sender, EventArgs e)
        {
            Console.WriteLine("dataGridView1[0, 0].Valueで値を呼び出せる");
            Console.WriteLine(dataGridView1[0, 0].Value);

            // XML file -> List<LensData>
            List<LensData> list = new List<LensData>();

            XmlSerializer serializer = new XmlSerializer(typeof(List<LensData>));
            using (var fs = new FileStream("lensdata.xml", FileMode.Open))
            {
                list =(List<LensData>)serializer.Deserialize(fs);
            }
            // List<LensData> -> DataGridView
            dataGridView1.RowCount = list.Count()+1; // いったん行数を伝える

            for (int row=0; row < dataGridView1.RowCount -1; row++)
            {
                LensData obj = list[row];

                dataGridView1[0, row].Value= obj.R;
                dataGridView1[1, row].Value = obj.D;
                dataGridView1[2, row].Value = obj.N;
                dataGridView1[3, row].Value = obj.Aperture;

            }

        }

        // Dataのセーブ
        private void button2_Click_1(object sender, EventArgs e)
        {
            // DATA Grid から値を吸い出して　LensDataにのリスト代入する
            // https://www.youtube.com/watch?v=QEOcrmcVCX0&list=WL&index=26

            List<LensData> list = new List<LensData>();

            for (int row = 0; row < dataGridView1.RowCount - 1; row++)
            {
                LensData obj = new LensData();
                obj.R = Convert.ToDouble(dataGridView1[0, row].Value);
                obj.D = Convert.ToDouble(dataGridView1[1, row].Value);
                obj.N = Convert.ToDouble(dataGridView1[2, row].Value);
                obj.Aperture = Convert.ToDouble(dataGridView1[3, row].Value);

                list.Add(obj);
            }
            //　データリストをXMLに吐き出す
            XmlSerializer serializer = new XmlSerializer(typeof(List<LensData>));
            using (var fs = new FileStream("lensdata.xml", FileMode.Create))
            {
                serializer.Serialize(fs, list);
            }

            
        }


        // お絵かき関数
        // https://dobon.net/vb/dotnet/graphics/drawline.html
        public void Line_Something()
        {
            // 画用紙みたいなもの？
            Bitmap canvas = new Bitmap(pictureBox2.Width, pictureBox2.Height);
            
            //　実際に筆やペンを設置して自動で動作する機械みたいなイメージ？
            //  ここでは、canvasに書きに行け！みたいな感じ？
            Graphics g = Graphics.FromImage(canvas);

            //　実際にペンで書く
            g.DrawLine(Pens.RoyalBlue, 10, 20, 50, 75);

            //　四角を書く
            g.DrawRectangle(Pens.Green, 20, 40, 50, 75);

            //　自動お絵かきマシーンのお片付け
            g.Dispose();

            //　canvasをPictureBox2に張り付ける
            pictureBox2.Image = canvas;

        }

        private void pictureBox2_Click(object sender, EventArgs e)
        {

        }

        private void SaveCSV_Click(object sender, EventArgs e)
        {
            // DataGridの内容をCSVファイルに出力するボタンを作る
            // https://qiita.com/nkimra/items/c8f42a9e2f3fda30414d

            //https://docs.microsoft.com/ja-jp/dotnet/api/microsoft.win32.savefiledialog?view=windowsdesktop-6.0
            // PresentationFramwork 追加する
            Microsoft.Win32.SaveFileDialog dlg = new Microsoft.Win32.SaveFileDialog();
            dlg.InitialDirectory = System.Environment.GetFolderPath(Environment.SpecialFolder.Personal);
            dlg.Title = "保存先のファイルを選択してください";
            dlg.Filter = "CSVファイル(*.csv)|*.csv";


            //　data = (ObservableCollection<User>)dataGridView1.DataSource;  //dataGridView1の値をdataに代入

            if (dlg.ShowDialog()==true) 
            {
                try
                {
                    using (var sw = new System.IO.StreamWriter(dlg.FileName, false, System.Text.Encoding.GetEncoding("Shift_JIS")))
                    {
                        //ダブルクォーテーションで囲む
                        Func<string, string> dqot = (str) => { return "\"" + str.Replace("\"", "\"\"") + "\""; };
                        foreach (User d in data)
                            sw.WriteLine(dqot(d.Name) + "," + dqot(d.Place));

                    }
                    System.Windows.Forms.MessageBox.Show("保存しました");
                }
                catch(SystemException ex)
                {
                    System.Console.WriteLine(ex.Message);
                }
                    
            }
        }

        private void OpenCSV_Click(object sender, EventArgs e)
        {
            // (参考サイト)[ VB.NET / C# ] CSVファイルからデータグリッドビューにデータをインポートする
            // https://hensa40.cutegirl.jp/archives/785
            // 注意点：MicrosoftBisualBasicの参照を追加する
            Microsoft.Win32.OpenFileDialog dlg = new Microsoft.Win32.OpenFileDialog();
            dlg.InitialDirectory = System.Environment.GetFolderPath(Environment.SpecialFolder.Personal);
            dlg.Title = "所望のファイルを選択してください";
            dlg.Filter = "CSVファイル(*.csv)|*.csv";

            dlg.ShowDialog();
            System.String str =dlg.FileName;
            // System.String str = @"C:\Users\13273_Yamazaki\source\repos\WindowsFormsApp1_CS\WindowsFormsApp1_CS\bin\x64\Release\lensdata.csv";

            if (String.IsNullOrEmpty(str))
            {
                System.Windows.Forms.MessageBox.Show("ファイルが読み込まれませんでした");
                return;
            }
            else { 

            TextFieldParser parser = new TextFieldParser(str, Encoding.GetEncoding("Shift_JIS"));
            parser.TextFieldType = FieldType.Delimited;
            parser.SetDelimiters(","); // 区切り文字はコンマ

            // データをすべてクリア
            dataGridView1.Rows.Clear();

            int i = 0;
                while (!parser.EndOfData)
                {
                    // 1行目だけは読み込まない（ヘッダーはフォーム初期化で記入）
                    if (i == 0)
                    {
                        string[] row = parser.ReadFields(); // 1行目だけ空読みさせる
                    }
                    else
                    {
                        string[] row = parser.ReadFields(); // 1行読み込み
                        // 読み込んだデータ(1行をDataGridViewに表示する)
                        dataGridView1.Rows.Add(row);
                    }
                    i += 1; // インクリメント
                }
            }
        }

        private void CSVsaveTest_Click(object sender, EventArgs e)
        {
            // ファイルダイアログで、ファイルパスを取得する
            Microsoft.Win32.SaveFileDialog dlg = new Microsoft.Win32.SaveFileDialog();
            dlg.InitialDirectory = System.Environment.GetFolderPath(Environment.SpecialFolder.Personal);
            dlg.Title = "保存先のファイルを選択してください";
            dlg.Filter = "CSVファイル(*.csv)|*.csv";
            string path = dlg.FileName;
            string strData = "";

            if (dlg.ShowDialog() == true)
            {
                try
                {
                    using (var sw = new System.IO.StreamWriter(dlg.FileName, false, System.Text.Encoding.GetEncoding("Shift_JIS")))
                    {
                        // 1行目（ヘッダー書き込み）
                        strData = "R" + ","
                            + "D" + ","
                            + "N" + ","
                            + "Aperture" + ",";
                        sw.WriteLine(strData);
                        // レンズデータ書き込み
                        for (int row = 0; row < dataGridView1.RowCount - 1; row++)
                        {
                            strData=Convert.ToDouble(dataGridView1[0, row].Value)+ ","
                            + Convert.ToDouble(dataGridView1[1, row].Value) + ","
                            + Convert.ToDouble(dataGridView1[2, row].Value) + ","
                            + Convert.ToDouble(dataGridView1[3, row].Value);

                            sw.WriteLine(strData);
                        }
                        // sw.Close();  //消す作業が必要？？？
                    }
                    // System.Windows.Forms.MessageBox.Show(path); // なぜ何も出てこない？
                    System.Windows.Forms.MessageBox.Show("保存しました");
                }
                catch (SystemException ex)
                {
                    System.Console.WriteLine(ex.Message);
                }
            }
        }

        private void Draw_Rect_Click(object sender, EventArgs e)
        {
            // パラメータの変更
            x_Rect = int.Parse(textBox_width.Text);
            y_Rect = int.Parse(textBox_height.Text);
            
            // pictureBox2.Paint += pictureBox2_Paint;
            pictureBox2.Refresh();
        }

        private void Draw_Rect_Paint(object sender, PaintEventArgs e)
        {
            

            
        }

        private void pictureBox2_Paint(object sender, PaintEventArgs e)
        {
            // ピクチャーボックスを描画するときの命令で e.Graphics が描画先の装置になる。
            var g = e.Graphics;

            // パラメータの設定
            // x_Rect = 111;
            // y_Rect = 35;

            g.DrawRectangle(
                new Pen(Color.FromArgb(255, 0, 0), 2),  // 赤色(RGBで指定)、線の太さ2
                10, 30, x_Rect, y_Rect // 長方形の左上と右下の座標を指定
                );
        }

        private void AS_Dist図_Click(object sender, EventArgs e)
        {

        }

        private void pictureBox1_Click(object sender, EventArgs e)
        {

        }

        private void pictureBox3_Click(object sender, EventArgs e)
        {

        }
    }

    // データをclassで定義する
    //　listにしたものをXML化する
    public class LensData
    {
        public double R;
        public double D;
        public double N;
        public double Aperture;
    }

    class User
    {
        public string Name { get; set; }
        public string Place { get; set; }

        public User(string name, string place)
        {
            this.Name = name;
            this.Place = place;
        }
    }
}
