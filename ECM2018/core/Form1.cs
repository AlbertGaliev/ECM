using ECM2018.core.edwards;
using ECM2018.core.montgomery;
using ECM2018.core.weierstrass;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Diagnostics;
using System.Threading;
using System.Threading.Tasks;
using System.Windows.Forms;
using static ECM2018.core.ECMWorker;
using ECM2018.core;
using ECM2018.core.algorithms;
using System.IO;

namespace ECM2018
{
    public partial class Form1 : Form
    {
        ECMParams param = new ECMParams();
        ECMWorker worker;
        Thread ecmThread;

        public Form1()
        {
            InitializeComponent();
            comboBox1.SelectedIndex = 0;
            comboBox2.SelectedIndex = 0;
            comboBox3.SelectedIndex = 0;
            comboBox4.SelectedIndex = 0;
            comboBox5.SelectedIndex = 0;
            button2.Enabled = false;
        }

        private void comboBox1_SelectedIndexChanged(object sender, EventArgs e)//сменили вид уравнения кривой
        {
            comboBox2.Items.Clear();
            comboBox2.Text = "";
            if (comboBox1.SelectedIndex == 0)
            {
                param.curveType = CurveType.WeierstrassCurve;
                comboBox2.Items.Add("Аффинные");
                comboBox2.Items.Add("Проективные");
                comboBox2.Items.Add("Якобиановы");
            }
            if (comboBox1.SelectedIndex == 1)
            {
                param.curveType = CurveType.MontgomeryCurve;
                comboBox2.Items.Add("Проективные XZ");
            }
            if (comboBox1.SelectedIndex == 2)
            {
                param.curveType = CurveType.TwistedEdwardsCurve;
                comboBox2.Items.Add("Проективные");
                comboBox2.Items.Add("Инвертированные");
                comboBox2.Items.Add("Расширенные");
            }
            comboBox4.Items.Clear();
            if (comboBox1.SelectedIndex == 0)
            {
                comboBox4.Items.Add("Случайные");
            }
            if (comboBox1.SelectedIndex == 1)
            {
                comboBox4.Items.Add("Случайные");
                comboBox4.Items.Add("Z/6Z");
            }
            if (comboBox1.SelectedIndex == 2)
            {
                comboBox4.Items.Add("Случайные");
                comboBox4.Items.Add("Z/12Z");
                comboBox4.Items.Add("Z/2Z x Z/8Z");
            }
            comboBox2.SelectedIndex = 0;
            comboBox4.SelectedIndex = 0;
            param.torsionType = TorsionType.RandomCurve;
        }

        private void comboBox3_SelectedIndexChanged(object sender, EventArgs e)//зафиксировали параметры В
        {
            if (comboBox3.SelectedIndex == 0)
            {
                param.autoApproximation = false;
                textBox1.Enabled = true;
                textBox2.Enabled = true;
            }
            if (comboBox3.SelectedIndex == 1)
            {
                param.autoApproximation = true;
                textBox1.Enabled = false;
                textBox2.Enabled = false;
            }
        }

        private void button1_Click(object sender, EventArgs e)//запустили алгоритм ЕСМ в отдельном потоке
        {
            ecmThread = new Thread(Lenstra);
            ecmThread.Start();
            new Thread(CountTime).Start(label7);
        }

        private void comboBox2_SelectedIndexChanged(object sender, EventArgs e)//сменили системы координат
        {
            switch (comboBox1.SelectedIndex)
            {
                case 0:
                    switch (comboBox2.SelectedIndex)
                    {
                        case 0:
                            param.coordinatesType = CoordinatesType.Affine;
                            break;
                        case 1:
                            param.coordinatesType = CoordinatesType.Projective;
                            break;
                        case 2:
                            param.coordinatesType = CoordinatesType.Jacobian;
                            break;
                    }
                    break;
                case 1:
                    switch (comboBox2.SelectedIndex)
                    {
                        case 0:
                            param.coordinatesType = CoordinatesType.Projective;
                            break;
                    }
                    break;
                case 2:
                    switch (comboBox2.SelectedIndex)
                    {
                        case 0:
                            param.coordinatesType = CoordinatesType.Projective;
                            break;
                        case 1:
                            param.coordinatesType = CoordinatesType.Inverted;
                            break;
                        case 2:
                            param.coordinatesType = CoordinatesType.Extended;
                            break;
                        case 3:
                            param.coordinatesType = CoordinatesType.Completed;
                            break;
                    }
                    break;
            }

        }

        private void comboBox4_SelectedIndexChanged(object sender, EventArgs e)//сменили кручение
        {
            switch (comboBox1.SelectedIndex)
            {
                case 0:
                    switch (comboBox4.SelectedIndex)
                    {
                        case 0:
                            param.torsionType = TorsionType.RandomCurve;
                            break;
                    }
                    break;
                case 1:
                    switch (comboBox4.SelectedIndex)
                    {
                        case 0:
                            param.torsionType = TorsionType.RandomCurve;
                            break;
                        case 1:
                            param.torsionType = TorsionType.Curve6;
                            break;
                    }
                    break;
                case 2:
                    switch (comboBox4.SelectedIndex)
                    {
                        case 0:
                            param.torsionType = TorsionType.RandomCurve;
                            break;
                        case 1:
                            param.torsionType = TorsionType.Curve12;
                            break;
                        case 2:
                            param.torsionType = TorsionType.Curve2x8;
                            break;
                    }
                    break;
            }

        }

        private void Lenstra()
        {
            BigInteger res, n;
            try
            {
                res = 0;
                textBox5.Text = "";
                label7.Text = "";
                button1.Enabled = false;
                button2.Enabled = true;
                n = BigInteger.Parse(textBox4.Text);
                param.B1 = Int32.Parse(textBox1.Text);
                param.B2 = Int32.Parse(textBox2.Text);
                param.curves = Int32.Parse(textBox3.Text);

            }
            catch (Exception)
            {
                MessageBox.Show("ВВЕДИТЕ КОРРЕКТНЫЕ ДАННЫЕ!");
                Thread.Sleep(1000);
                button1.Enabled = true;
                button2.Enabled = false;
                label7.Text = "";
                return;
            }
            worker = new ECMWorker(param);
            var sw = new Stopwatch();
            sw.Start();
            if (comboBox5.SelectedIndex == 0)
            {
                var result = worker.GetFactor(n);
                res = result;
                if (result == 0)
                    textBox5.Text = "Не удалось найти делитель";
                else
                    textBox5.Text = result.ToString();
            }
            else
            {
                var result = worker.GetAllFactors(n);
                string str = "";
                if (result.Length == 0)
                    textBox5.Text = "Не удалось найти делители";
                else
                {
                    for (int i = 0; i < result.Length; i++)
                    {
                        str += result[i].ToString();
                        str += "*";
                    }
                    str = str.Remove(str.Length - 1, 1);
                    textBox5.Text = str.ToString();
                }
            }
            sw.Stop();
            label7.Text = sw.Elapsed.TotalMilliseconds.ToString() + " мс";
            button1.Enabled = true;
            button2.Enabled = false;
        }

        private void button2_Click(object sender, EventArgs e)// остановили алгоритм
        {
            button2.Enabled = false;
            worker.StopFactoring();
        }

        private void button3_Click(object sender, EventArgs e)//генерация РСА
        {
            textBox7.Text = "";
            textBox8.Text = "";
            textBox9.Text = "";
            label8.Text = "";
            label9.Text = "";
            label10.Text = "";
            int length;
            try
            {
                 length = Int32.Parse(textBox6.Text);
            }
            catch (Exception)
            {
                MessageBox.Show("ВВЕДИТЕ КОРРЕКТНЫЕ ДАННЫЕ!");
                return;
            }
            BigInteger p = Algorithms.PrimeGenerator(length);
            BigInteger q = Algorithms.PrimeGenerator(length);
            BigInteger N = p * q;
            textBox7.Text = p.ToString();
            textBox8.Text = q.ToString();
            textBox9.Text = N.ToString();
            label8.Text = p.ToString().Length.ToString() + "-значное число";
            label9.Text = q.ToString().Length.ToString() + "-значное число";
            label10.Text = N.ToString().Length.ToString() + "-значное число";
        }

        private void button4_Click(object sender, EventArgs e)//проверка простоты
        {
            try
            {
                label17.Text = "";
                if (Algorithms.MillerRabinTest(BigInteger.Parse(textBox10.Text)))
                {
                    label17.Text = "простое.";
                }
                else
                {
                    label17.Text = "составное.";
                }

            }
            catch (Exception)
            {
                MessageBox.Show("ВВЕДИТЕ КОРРЕКТНЫЕ ДАННЫЕ!");
                return;
            }
        }

        private void CountTime(object labelObject)
        {
            var label = labelObject as Label;
            int time = 0;
            while (ecmThread.IsAlive)
            {
                label.Text = time + " с";
                Thread.Sleep(1000);
                time++;
            }
        }

    }
}
