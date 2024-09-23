using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms; //Net framework 4.8

namespace MetNumericosC
{
    internal class MetSolEcuaciones
    {
        public int NumMaxIter;
        public float ErrorMaximo;
        float AproxRaiz;
        Complex AproxRaiC;

        public bool MetBiseccion(float a, float b, ref DataGridView dgvResultado)
        {
            float c, ErrorAct;
            int i;

            if (Func(a) * Func(b) > 0)
            {
                MessageBox.Show("No se puede aplicar el método de bisección", "ERROR", MessageBoxButtons.OK,
                    MessageBoxIcon.Error);
                return false;
            }

            dgvResultado.Rows.Clear();

            dgvResultado.Columns.Add("Iteracion", "i");
            dgvResultado.Columns.Add("Valor_a", "a");
            dgvResultado.Columns.Add("Valor_c", "c");
            dgvResultado.Columns.Add("Valor_b", "b");
            dgvResultado.Columns.Add("F_a", "f(a)");
            dgvResultado.Columns.Add("F_c", "c");
            dgvResultado.Columns.Add("F_b", "b");
            dgvResultado.Columns.Add("Error", "error");



            i = 1; //falta ver si este VA
            while (i <= NumMaxIter)
            {
                c = (a + b) / 2;
                ErrorAct = (b - a) / 2; //error actual

                //llenando las filas
                dgvResultado.Rows.Add();
                dgvResultado.Rows[i - 1].Cells[0].Value = i;
                dgvResultado.Rows[i - 1].Cells[1].Value = a;
                dgvResultado.Rows[i - 1].Cells[2].Value = c;
                dgvResultado.Rows[i - 1].Cells[3].Value = b;
                dgvResultado.Rows[i - 1].Cells[4].Value = Func(a) > 0 ? "+" : "-";
                dgvResultado.Rows[i - 1].Cells[5].Value = Func(c) > 0 ? "+" : "-";
                dgvResultado.Rows[i - 1].Cells[6].Value = Func(b) > 0 ? "+" : "-";
                dgvResultado.Rows[i - 1].Cells[7].Value = ErrorAct;


                if (ErrorAct <= ErrorMaximo)
                {
                    MessageBox.Show("Se obtuvo la raíz con el error deseado. Raíz=" + c.ToString(), "AVISO", MessageBoxButtons.OK,
                    MessageBoxIcon.Information);
                    return true;
                }
                if (Func(a) * Func(c) < 0)
                    b = c;
                else
                    a = c;

                i++;
            }
            MessageBox.Show("No se pudo obtener la aprox con error deseado.", "AVISO", MessageBoxButtons.OK,
                    MessageBoxIcon.Warning);
            return false;
        } // se termina método bisección 


        private float Func(float x)
        {
            float r;
            r = (float)(Math.Pow((double)x, (double)2.0) - 3.0);
            return r;
        }
        dgvResultado.Rows.Clear(); 

            dgvResultado.Columns.Add("iteracion", "i");
            dgvResultado.Columns.Add("Valor_a", "a");
            dgvResultado.Columns.Add("Valor_b", "b");
            dgvResultado.Columns.Add("Valor_c", "c");
            dgvResultado.Columns.Add("F_a", "f(a)");
            dgvResultado.Columns.Add("F_c", "f(c)");
            dgvResultado.Columns.Add("F_b", "f(b)");
            dgvResultado.Columns.Add("Error", "error");


            i = 1;
            while (i <= NumMaxIter)
            {
                // p = a - f(a)+(b-a)/(f(b)-f(a)); 
                c = a-Func(a)+(b-a)/Func()
                ErrorAct = (b-a) / 2;

                dgvResultado.Rows.Add(); 
                dgvResultado.Rows[i - 1].Cells[2].Value = b; 
                
            }
    /* 
        c= f(p2)
        b= (p0-p2)^2 * (f(p1) - f(p2)) - (p1 - p2)^2 * (f(p0)-f(p1))
        b = b / ((p0-p2)*(p1-p2)*(p0-p1))


     */

    public bool MetMuller(Complex p0, Complex p1, Complex p2, ref DataGridView dgvResultado)
    {
        float ErrorAct;
        Complex p = new Complex(), a = new Complex(), b = new Complex(), c = new Complex();
        Complex Fp1MFp2 = new Complex();
        Complex Fp0MFp2 = new Complex();
        Complex p0Mp2 = new Complex();
        Complex p0Mp1 = new Complex();
        Complex p1Mp2 = new Complex();
        Complex divisor = new Complex();
        Complex D = new Complex();
        Complex E = new Complex();
        Complex Aux = new Complex();
        int i;


        dgvResultado.Rows.Clear();
        dgvResultado.Columns.Clear();

        dgvResultado.Columns.Add("iteracion", "i");
        dgvResultado.Columns.Add("Valor_p0", "p0");
        dgvResultado.Columns.Add("Valor_p1", "p1");
        dgvResultado.Columns.Add("Valor_p2", "p2");
        dgvResultado.Columns.Add("p", "p");
        dgvResultado.Columns.Add("Error", "error");

        i = 1;
        while (i <= NumMaxIter)
        {
            p0Mp2 = Complex.Subtract(p0, p2);
            p0Mp1 = Complex.Subtract(p0, p1);
            p1Mp2 = Complex.Subtract(p1, p2);
            Fp1MFp2 = Complex.Subtract(FuncCompleja(p1), FuncCompleja(p2));
            Fp0MFp2 = Complex.Subtract(FuncCompleja(p0), FuncCompleja(p2));

            divisor = Complex.Multiply(p0Mp2, p1Mp2);
            divisor = Complex.Multiply(divisor, p0Mp1);

            c = FuncCompleja(p2);

            b = Complex.Subtract(Complex.Multiply(Complex.Pow(p0Mp2, 2), Fp1MFp2),
               Complex.Multiply(Complex.Pow(p1Mp2, 2), Fp0MFp2));
            b = Complex.Divide(b, divisor);

            a = Complex.Subtract(Complex.Multiply(p1Mp2, Fp0MFp2), Complex.Multiply(p0Mp2, Fp1MFp2));
            a = Complex.Divide(a, divisor);

            D = Complex.Sqrt(Complex.Subtract(Complex.Pow(b, 2), Complex.Multiply(new Complex(4, 0), Complex.Multiply(a, c))));
            if (Complex.Abs(Complex.Subtract(b, D)) < Complex.Abs(Complex.Add(b, D)))
                E = Complex.Add(b, D);
            p = Complex.Subtract(p2, Complex.Divide(Complex.Multiply(new Complex(2, 0), c), E));

            ErrorAct = (float)Complex.Abs(Complex.Divide(Complex.Multiply(new Complex(2, 0), c), E));

            dgvResultado.Rows.Add();
            dgvResultado.Rows[i - 1].Cells[0].Value = i;
            dgvResultado.Rows[i - 1].Cells[1].Value = p0;
            dgvResultado.Rows[i - 1].Cells[2].Value = p1;
            dgvResultado.Rows[i - 1].Cells[3].Value = p2;
            dgvResultado.Rows[i - 1].Cells[4].Value = p;

            dgvResultado.Rows[i - 1].Cells[5].Value = ErrorAct;

            if (ErrorAct <= ErrorMaximo)
            {
                MessageBox.Show("Se obtuvo la raíz con el error deseado. Raiz" + p.ToString(), "AVISO", MessageBoxButtons.OK,
                    MessageBoxIcon.Information);
                return true;
            }
            p0 = p1;
            p1 = p2;
            p2 = p;
            i++;
        }

        MessageBox.Show("No se pudo obtener la aprox con error deseado.", "AVISO", MessageBoxButtons.OK,
                    MessageBoxIcon.Warning);
        return false;
    }

    private float Func(float x)
    {
        float r;
        r = (float)(Math.Pow((double)x, (double)2.0 - 3.0);
        return r;
    }

    private Complex FuncCompleja(Complex x)
    {
        Complex r = new Complex();
        r = Complex.Pow(x, 2) - new Complex(3, 0);
        return r;
    }
}  
}
