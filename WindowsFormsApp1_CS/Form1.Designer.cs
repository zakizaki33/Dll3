namespace WindowsFormsApp1_CS
{
    partial class Form1
    {
        /// <summary>
        /// 必要なデザイナー変数です。
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// 使用中のリソースをすべてクリーンアップします。
        /// </summary>
        /// <param name="disposing">マネージド リソースを破棄する場合は true を指定し、その他の場合は false を指定します。</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows フォーム デザイナーで生成されたコード

        /// <summary>
        /// デザイナー サポートに必要なメソッドです。このメソッドの内容を
        /// コード エディターで変更しないでください。
        /// </summary>
        private void InitializeComponent()
        {
            this.button1 = new System.Windows.Forms.Button();
            this.dataSet1 = new System.Data.DataSet();
            this.label1 = new System.Windows.Forms.Label();
            this.dataGridView1 = new System.Windows.Forms.DataGridView();
            this.textBox1 = new System.Windows.Forms.TextBox();
            this.tabControl1 = new System.Windows.Forms.TabControl();
            this.LensView = new System.Windows.Forms.TabPage();
            this.SA図 = new System.Windows.Forms.TabPage();
            this.pictureBox1 = new System.Windows.Forms.PictureBox();
            this.AS_Dist図 = new System.Windows.Forms.TabPage();
            this.textBox_width = new System.Windows.Forms.TextBox();
            this.label3 = new System.Windows.Forms.Label();
            this.label2 = new System.Windows.Forms.Label();
            this.textBox_height = new System.Windows.Forms.TextBox();
            this.Draw_Rect = new System.Windows.Forms.Button();
            this.pictureBox2 = new System.Windows.Forms.PictureBox();
            this.Open = new System.Windows.Forms.Button();
            this.button2 = new System.Windows.Forms.Button();
            this.SaveCSV = new System.Windows.Forms.Button();
            this.OpenCSV = new System.Windows.Forms.Button();
            this.CSVsaveTest = new System.Windows.Forms.Button();
            ((System.ComponentModel.ISupportInitialize)(this.dataSet1)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.dataGridView1)).BeginInit();
            this.tabControl1.SuspendLayout();
            this.LensView.SuspendLayout();
            this.SA図.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.pictureBox1)).BeginInit();
            this.AS_Dist図.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.pictureBox2)).BeginInit();
            this.SuspendLayout();
            // 
            // button1
            // 
            this.button1.Location = new System.Drawing.Point(545, 209);
            this.button1.Name = "button1";
            this.button1.Size = new System.Drawing.Size(66, 43);
            this.button1.TabIndex = 0;
            this.button1.Text = "Calc";
            this.button1.UseVisualStyleBackColor = true;
            this.button1.Click += new System.EventHandler(this.button1_Click);
            // 
            // dataSet1
            // 
            this.dataSet1.DataSetName = "NewDataSet";
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(544, 297);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(110, 12);
            this.label1.TabIndex = 1;
            this.label1.Text = "Focal Length = ?????";
            this.label1.Click += new System.EventHandler(this.label1_Click);
            // 
            // dataGridView1
            // 
            this.dataGridView1.AllowDrop = true;
            this.dataGridView1.ColumnHeadersHeightSizeMode = System.Windows.Forms.DataGridViewColumnHeadersHeightSizeMode.AutoSize;
            this.dataGridView1.Location = new System.Drawing.Point(52, 53);
            this.dataGridView1.Name = "dataGridView1";
            this.dataGridView1.RowTemplate.Height = 21;
            this.dataGridView1.Size = new System.Drawing.Size(442, 272);
            this.dataGridView1.TabIndex = 2;
            this.dataGridView1.DragDrop += new System.Windows.Forms.DragEventHandler(this.dataGridView1_DragDrop);
            this.dataGridView1.DragEnter += new System.Windows.Forms.DragEventHandler(this.dataGridView1_DragEnter);
            // 
            // textBox1
            // 
            this.textBox1.Font = new System.Drawing.Font("メイリオ", 14.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(128)));
            this.textBox1.Location = new System.Drawing.Point(32, 37);
            this.textBox1.Name = "textBox1";
            this.textBox1.Size = new System.Drawing.Size(195, 36);
            this.textBox1.TabIndex = 3;
            this.textBox1.TextChanged += new System.EventHandler(this.textBox1_TextChanged);
            // 
            // tabControl1
            // 
            this.tabControl1.Controls.Add(this.LensView);
            this.tabControl1.Controls.Add(this.SA図);
            this.tabControl1.Controls.Add(this.AS_Dist図);
            this.tabControl1.Location = new System.Drawing.Point(52, 368);
            this.tabControl1.Name = "tabControl1";
            this.tabControl1.SelectedIndex = 0;
            this.tabControl1.Size = new System.Drawing.Size(785, 224);
            this.tabControl1.TabIndex = 4;
            // 
            // LensView
            // 
            this.LensView.Controls.Add(this.textBox1);
            this.LensView.Location = new System.Drawing.Point(4, 22);
            this.LensView.Name = "LensView";
            this.LensView.Padding = new System.Windows.Forms.Padding(3);
            this.LensView.Size = new System.Drawing.Size(777, 198);
            this.LensView.TabIndex = 0;
            this.LensView.Text = "LensView";
            this.LensView.UseVisualStyleBackColor = true;
            // 
            // SA図
            // 
            this.SA図.Controls.Add(this.pictureBox1);
            this.SA図.Location = new System.Drawing.Point(4, 22);
            this.SA図.Name = "SA図";
            this.SA図.Padding = new System.Windows.Forms.Padding(3);
            this.SA図.Size = new System.Drawing.Size(777, 198);
            this.SA図.TabIndex = 1;
            this.SA図.Text = "SA図";
            this.SA図.UseVisualStyleBackColor = true;
            // 
            // pictureBox1
            // 
            this.pictureBox1.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle;
            this.pictureBox1.Location = new System.Drawing.Point(38, 18);
            this.pictureBox1.Name = "pictureBox1";
            this.pictureBox1.Size = new System.Drawing.Size(161, 162);
            this.pictureBox1.SizeMode = System.Windows.Forms.PictureBoxSizeMode.Zoom;
            this.pictureBox1.TabIndex = 0;
            this.pictureBox1.TabStop = false;
            this.pictureBox1.DragDrop += new System.Windows.Forms.DragEventHandler(this.pictureBox1_DragDrop);
            this.pictureBox1.DragEnter += new System.Windows.Forms.DragEventHandler(this.pictureBox1_DragEnter);
            // 
            // AS_Dist図
            // 
            this.AS_Dist図.Controls.Add(this.textBox_width);
            this.AS_Dist図.Controls.Add(this.label3);
            this.AS_Dist図.Controls.Add(this.label2);
            this.AS_Dist図.Controls.Add(this.textBox_height);
            this.AS_Dist図.Controls.Add(this.Draw_Rect);
            this.AS_Dist図.Controls.Add(this.pictureBox2);
            this.AS_Dist図.Location = new System.Drawing.Point(4, 22);
            this.AS_Dist図.Name = "AS_Dist図";
            this.AS_Dist図.Padding = new System.Windows.Forms.Padding(3);
            this.AS_Dist図.Size = new System.Drawing.Size(777, 198);
            this.AS_Dist図.TabIndex = 2;
            this.AS_Dist図.Text = "AS_DIST図";
            this.AS_Dist図.UseVisualStyleBackColor = true;
            this.AS_Dist図.Click += new System.EventHandler(this.AS_Dist図_Click);
            // 
            // textBox_width
            // 
            this.textBox_width.Location = new System.Drawing.Point(576, 126);
            this.textBox_width.Name = "textBox_width";
            this.textBox_width.Size = new System.Drawing.Size(100, 19);
            this.textBox_width.TabIndex = 5;
            // 
            // label3
            // 
            this.label3.AutoSize = true;
            this.label3.Location = new System.Drawing.Point(576, 94);
            this.label3.Name = "label3";
            this.label3.Size = new System.Drawing.Size(17, 12);
            this.label3.TabIndex = 4;
            this.label3.Text = "幅";
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Location = new System.Drawing.Point(574, 21);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(25, 12);
            this.label2.TabIndex = 3;
            this.label2.Text = "高さ";
            // 
            // textBox_height
            // 
            this.textBox_height.Location = new System.Drawing.Point(574, 55);
            this.textBox_height.Name = "textBox_height";
            this.textBox_height.Size = new System.Drawing.Size(100, 19);
            this.textBox_height.TabIndex = 2;
            // 
            // Draw_Rect
            // 
            this.Draw_Rect.Location = new System.Drawing.Point(595, 158);
            this.Draw_Rect.Name = "Draw_Rect";
            this.Draw_Rect.Size = new System.Drawing.Size(79, 20);
            this.Draw_Rect.TabIndex = 1;
            this.Draw_Rect.Text = "矩形を描く";
            this.Draw_Rect.UseVisualStyleBackColor = true;
            this.Draw_Rect.Click += new System.EventHandler(this.Draw_Rect_Click);
            this.Draw_Rect.Paint += new System.Windows.Forms.PaintEventHandler(this.Draw_Rect_Paint);
            // 
            // pictureBox2
            // 
            this.pictureBox2.Location = new System.Drawing.Point(31, 21);
            this.pictureBox2.Name = "pictureBox2";
            this.pictureBox2.Size = new System.Drawing.Size(486, 157);
            this.pictureBox2.TabIndex = 0;
            this.pictureBox2.TabStop = false;
            this.pictureBox2.Click += new System.EventHandler(this.pictureBox2_Click);
            this.pictureBox2.Paint += new System.Windows.Forms.PaintEventHandler(this.pictureBox2_Paint);
            // 
            // Open
            // 
            this.Open.Location = new System.Drawing.Point(550, 61);
            this.Open.Name = "Open";
            this.Open.Size = new System.Drawing.Size(61, 31);
            this.Open.TabIndex = 5;
            this.Open.Text = "Open";
            this.Open.UseVisualStyleBackColor = true;
            this.Open.Click += new System.EventHandler(this.button2_Click);
            // 
            // button2
            // 
            this.button2.Location = new System.Drawing.Point(550, 120);
            this.button2.Name = "button2";
            this.button2.Size = new System.Drawing.Size(61, 31);
            this.button2.TabIndex = 6;
            this.button2.Text = "Save";
            this.button2.UseVisualStyleBackColor = true;
            this.button2.Click += new System.EventHandler(this.button2_Click_1);
            // 
            // SaveCSV
            // 
            this.SaveCSV.Location = new System.Drawing.Point(707, 124);
            this.SaveCSV.Name = "SaveCSV";
            this.SaveCSV.Size = new System.Drawing.Size(96, 27);
            this.SaveCSV.TabIndex = 7;
            this.SaveCSV.Text = "save_CSV";
            this.SaveCSV.UseVisualStyleBackColor = true;
            this.SaveCSV.Click += new System.EventHandler(this.SaveCSV_Click);
            // 
            // OpenCSV
            // 
            this.OpenCSV.Location = new System.Drawing.Point(710, 61);
            this.OpenCSV.Name = "OpenCSV";
            this.OpenCSV.Size = new System.Drawing.Size(93, 24);
            this.OpenCSV.TabIndex = 8;
            this.OpenCSV.Text = "open_CSV";
            this.OpenCSV.UseVisualStyleBackColor = true;
            this.OpenCSV.Click += new System.EventHandler(this.OpenCSV_Click);
            // 
            // CSVsaveTest
            // 
            this.CSVsaveTest.Location = new System.Drawing.Point(711, 214);
            this.CSVsaveTest.Name = "CSVsaveTest";
            this.CSVsaveTest.Size = new System.Drawing.Size(100, 37);
            this.CSVsaveTest.TabIndex = 9;
            this.CSVsaveTest.Text = "saveCSVtest";
            this.CSVsaveTest.UseVisualStyleBackColor = true;
            this.CSVsaveTest.Click += new System.EventHandler(this.CSVsaveTest_Click);
            // 
            // Form1
            // 
            this.AllowDrop = true;
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 12F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(893, 604);
            this.Controls.Add(this.CSVsaveTest);
            this.Controls.Add(this.OpenCSV);
            this.Controls.Add(this.SaveCSV);
            this.Controls.Add(this.button2);
            this.Controls.Add(this.Open);
            this.Controls.Add(this.tabControl1);
            this.Controls.Add(this.dataGridView1);
            this.Controls.Add(this.label1);
            this.Controls.Add(this.button1);
            this.Name = "Form1";
            this.Text = "Form1";
            this.Load += new System.EventHandler(this.Form1_Load);
            ((System.ComponentModel.ISupportInitialize)(this.dataSet1)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.dataGridView1)).EndInit();
            this.tabControl1.ResumeLayout(false);
            this.LensView.ResumeLayout(false);
            this.LensView.PerformLayout();
            this.SA図.ResumeLayout(false);
            ((System.ComponentModel.ISupportInitialize)(this.pictureBox1)).EndInit();
            this.AS_Dist図.ResumeLayout(false);
            this.AS_Dist図.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)(this.pictureBox2)).EndInit();
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.Button button1;
        private System.Data.DataSet dataSet1;
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.DataGridView dataGridView1;
        private System.Windows.Forms.TextBox textBox1;
        private System.Windows.Forms.TabControl tabControl1;
        private System.Windows.Forms.TabPage LensView;
        private System.Windows.Forms.TabPage SA図;
        private System.Windows.Forms.TabPage AS_Dist図;
        private System.Windows.Forms.PictureBox pictureBox1;
        private System.Windows.Forms.Button Open;
        private System.Windows.Forms.Button button2;
        private System.Windows.Forms.PictureBox pictureBox2;
        private System.Windows.Forms.Button SaveCSV;
        private System.Windows.Forms.Button OpenCSV;
        private System.Windows.Forms.Button CSVsaveTest;
        private System.Windows.Forms.TextBox textBox_width;
        private System.Windows.Forms.Label label3;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.TextBox textBox_height;
        private System.Windows.Forms.Button Draw_Rect;
    }
}

