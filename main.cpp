#include <iostream>
#include "apmatrix.h"
#include <cmath>
#include <cstdlib>
#include <ctime>
using namespace std;
apmatrix<float> A,B,C,D,E,F,work;
apmatrix<float> c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10;
int counter = 0;
void EnterMatrix();
void PrintMatrix(apmatrix<float> selection);
apmatrix<float> TransposeMatrix(apmatrix<float> selection);
apmatrix<float> AddMatricies(apmatrix<float> selection1, apmatrix<float> selection2);
apmatrix<float> SubtractMatricies(apmatrix<float> selection1, apmatrix<float> selection2);
apmatrix<float> MultiplyMatricies(apmatrix<float> selection1, apmatrix<float> selection2);
apmatrix<float> RowReduceMatrix(apmatrix<float> in);
apmatrix<float> Invert(apmatrix<float> in);
apmatrix<float> createIdentity(int size);
apmatrix<float> augmentIdentity(apmatrix<float> in);
apmatrix<float> OLS(apmatrix<float> m, apmatrix<float> y);
apmatrix<float> ScalarMultiplyMatrix(apmatrix<float> in, float scalar);
apmatrix<float> VectorSubtract(apmatrix<float> in1, apmatrix<float> in2);
apmatrix<float> MakeQMatrix(apmatrix<float> a);
apmatrix<float> MakeRMatrix(apmatrix<float> a, apmatrix<float> Q);
apmatrix<float> QR(apmatrix<float> A, int times);
apmatrix<float> geteigenvalues(apmatrix<float> A, int times);
apmatrix<float> projection(apmatrix<float> thi, apmatrix<float> ontothat);
float dotprodvec(apmatrix<float> u, apmatrix<float> v);
float vectormagnitude(apmatrix<float> vec);
float pointpointDist(apmatrix<float> A, apmatrix<float> B);
apmatrix<float> pointlineDist(apmatrix<float> point, apmatrix<float> linePoint, apmatrix<float> directionVec);
void Swap(apmatrix<float> &in, int r1, int r2);
void LinearComb(apmatrix<float> &in, int r1, int lead, int leadcol);
void RowMultiply(apmatrix<float> &in, int lead, int leadcol);
void SetMatrix(apmatrix<float> target, apmatrix<float> newMatrix);
void testTheory(apmatrix<float> mat);
apmatrix<float> LoadMatrix(int matrixNumber);
float floatRand();
int main()
{
    int choice, cont;//, matrix_number;
    bool keepGoing = true;
    while(keepGoing)
    {
        cout<<"Welcome to Rohan's Matrix Machine \n\n";
        cout<<"What would you like to do? \n\n";
        cout<<"1. Enter a Matrix\n";
        cout<<"2. Transpose a Matrix\n";
        cout<<"3. Add Matricies\n";
        cout<<"4. Multiply Matricies\n";
        cout<<"5. Row Reduce a Matrix\n";
        cout<<"6. Invert a Matrix\n";
        cout<<"7. OLS Regression\n";
        cout<<"8. QR Algorithm\n";
        cin >> choice;
        cout<< "You picked option " <<choice <<"\n";
        if (choice == 1)
        {
            EnterMatrix();
        }
        else if (choice == 2)
        {
            int matrixNumber;
            cout<<"Which Matrix would you like to transpose? (1-6) \n";
            cin>>matrixNumber;
            work = LoadMatrix(matrixNumber);
            work = TransposeMatrix(work);
            PrintMatrix(work);
        }
        else if (choice == 3)
        {
            int first,second;
            cout<<"Pick your first matrix to add (1-6): \n";
            cin>>first;
            cout<<"pick your second matrix to add (1-6): \n";
            cin>>second;
            work = AddMatricies(LoadMatrix(first), LoadMatrix(second));
            PrintMatrix(work);
        }
        else if (choice == 4)
        {
            int first,second;
            cout<<"Pick your first matrix to multiply (1-6): \n";
            cin>>first;
            cout<<"pick your second matrix to multiply (1-6): \n";
            cin>>second;
            work = MultiplyMatricies(LoadMatrix(first), LoadMatrix(second));
            PrintMatrix(work);
        }
        else if (choice == 5)
        {
            int first;
            cout<<"Pick your matrix to row reduce (1-6): \n";
            cin>>first;
            apmatrix<float> work = LoadMatrix(first);
            work = RowReduceMatrix(work);
            PrintMatrix(work);
        }
        else if (choice == 6)
        {
            int first;
            cout<<"Pick your matrix to invert (must be square) (1-6): \n";
            cin>>first;
            apmatrix<float> work = LoadMatrix(first);
            work = Invert(work);
            PrintMatrix(work);
        }
        else if (choice == 7)
        {
            int size, order;
            cout<<"What order polynomial do you want your best fit to be? : \n";
            cin>>order;
            cout<<"How many data sets will you be adding? : \n";
            cin>>size;
            apmatrix<float> m,y;
            m.resize(size, order+1);
            y.resize(size, 1);
            for (int i = 0; i < size; i++)
            {
                float xin,yin;
                cout<<"Input your (x,y) pair: \n";
                cin>> xin;
                cin>> yin;
                y[i][0] = yin;
                for (int j = 0; j < order+1; j++)
                m[i][j] = pow(xin,j);
            }
            PrintMatrix(m);
            PrintMatrix(y);
            work = OLS(m,y);
            PrintMatrix(work);
        }
        else if (choice == 8)
        {
            int first, times;
            cout<<"Pick your matrix to find the eigenvalues of (1-6): \n";
            cin>>first;
            cout<<"How many times should the QR algorithm loop: \n";
            cin>>times;
            apmatrix<float> work = LoadMatrix(first);
            work = geteigenvalues(work, times);
            PrintMatrix(work);
        }
        else if (choice == 9)
        {
            PrintMatrix(A);
            //MakeQMatrix(A);
//            PrintMatrix(A);
//            PrintMatrix(B);
//            PrintMatrix(C);
//            cout<<"\n";
//            PrintMatrix(pointlineDist(A,B,C));
        }
        cout<<"Would you like to continue? 0 for no, 1 for yes: \n";
        cin>>cont;
        if (cont == 0)
        {
            keepGoing = false;
        }
        else
        {
            keepGoing = true;
        }
    }
}

void PrintMatrix(apmatrix<float> selection)
{
    for (int r = 0; r < selection.numrows(); r++)
    {
        for (int c = 0; c < selection.numcols(); c++)
        {
            cout<<selection[r][c]<<"\t";
        }
        cout<<"\n";
    }
}

void EnterMatrix()
{
    int row, col;
    apmatrix<float> matrix1;
    
    counter++;
    
    cout<<"How many rows? : ";
    cin>>row;
    cout<<"How many columns? : ";
    cin>>col;
    
    matrix1.resize(row, col);
    
    for(int r = 0; r < row; r++)
    {
        for(int c = 0; c < col; c++)
        {
            cout<<"Enter your row "<<r+1<<" column "<<c+1<<" entry: ";
            cin>>matrix1[r][c];
        }
    }
    if(counter == 1)
    {
        A = matrix1;
    }
    else if(counter == 2)
    {
        B = matrix1;
    }
    else if(counter == 3)
    {
        C = matrix1;
    }
    else if(counter == 4)
    {
        D = matrix1;
    }
    else if(counter == 5)
    {
        E = matrix1;
    }
    else if(counter == 6)
    {
        F = matrix1;
    }
    
    PrintMatrix(matrix1);
}

apmatrix<float> TransposeMatrix(apmatrix<float> selection)
{
    apmatrix<float> selectionTranspose;
    selectionTranspose.resize(selection.numcols(), selection.numrows());
    for (int i = 0; i < selection.numcols(); i++)
    {
        for (int j = 0; j < selection.numrows(); j++)
        {
            selectionTranspose[i][j] = selection[j][i];
        }
    }
    return selectionTranspose;
}

apmatrix<float> AddMatricies(apmatrix<float> selection1, apmatrix<float> selection2)
{
    apmatrix<float> sum;
    if (!(selection1.numrows() == selection2.numrows()) && (selection1.numcols() == selection2.numcols()))
    {
        cout<<"Invalid dimensions, could not add. Returning empty matrix";
    }
    else
    {
        sum.resize(selection1.numrows(), selection1.numcols());
        for (int i = 0; i < selection1.numrows(); i++)
        {
            for (int j = 0; j < selection1.numcols(); j++)
            {
                sum[i][j] = selection1[i][j] + selection2[i][j];
            }
        }
    }
    return sum;
}
apmatrix<float> SubtractMatricies(apmatrix<float> selection1, apmatrix<float> selection2)
{
    apmatrix<float> sum;
    if (!(selection1.numrows() == selection2.numrows()) && (selection1.numcols() == selection2.numcols()))
    {
        cout<<"Invalid dimensions, could not add. Returning empty matrix";
    }
    else
    {
        sum.resize(selection1.numrows(), selection1.numcols());
        for (int i = 0; i < selection1.numrows(); i++)
        {
            for (int j = 0; j < selection1.numcols(); j++)
            {
                sum[i][j] = selection1[i][j] - selection2[i][j];
            }
        }
    }
    return sum;
}
apmatrix<float> MultiplyMatricies(apmatrix<float> selection1, apmatrix<float> selection2)
{
    apmatrix<float> prod;
    if (!(selection1.numcols() == selection2.numrows()))
    {
        cout<<"Invalid dimensions, could not multiply. Returning empty matrix";
    }
    else
    {
        prod.resize(selection1.numrows(), selection2.numcols());
        for (int i = 0; i < prod.numrows(); i++)
        {
            for (int j = 0; j < prod.numcols(); j++)
            {
                for (int k = 0; k < selection1.numcols(); k++)
                {
                    prod[i][j] += selection1[i][k] * selection2[k][j];
                }
            }
        }
    }
    return prod;
}

apmatrix<float> LoadMatrix(int matrixNumber)
{
    apmatrix<float> Temp;
    
     if (matrixNumber == 1)
     {
     Temp = A;
     }
     else if (matrixNumber == 2)
     {
     Temp = B;
     }
     else if (matrixNumber == 3)
     {
         Temp = C;
     }
     else if (matrixNumber == 4)
     {
     Temp = D;
     }
     else if (matrixNumber == 5)
     {
     Temp = E;
     }
     else if (matrixNumber == 6)
     {
     Temp = F;
     }
    return Temp;
}

void SetMatrix(apmatrix<float> target, apmatrix<float> newMatrix)
{
    target = newMatrix;
}

void testTheory(apmatrix<float> mat)
{
    //apmatrix<float> tempq = MakeQMatrix(A);
    //apmatrix<float> tempr = MakeRMatrix(A, tempq);
    PrintMatrix(QR(A, 5));
}
apmatrix<float> RowReduceMatrix(apmatrix<float> in)
{
    apmatrix<float> red = in;
    int rows = red.numrows();
    int lead = 0;
    int leadcol = 0;
    while (lead < rows)
    {
        if(red[lead][leadcol] == 0)
        {
            for (int a = 1; a < rows - lead; a++)
            {
                if (lead + a < rows && red[lead+a][leadcol] != 0)
                {
                    Swap(red, lead, lead+a);
                    RowMultiply(red, lead, leadcol);
                    for (int i = 0; i < rows; i++)
                    {
                        if (i != lead)
                        {
                            LinearComb(red, i, lead, lead);
                        }
                    }
                }
                else if (lead+a >= rows)
                {
                    break;
                }
            }
            
        }
        if (red[lead][leadcol] == 1)
        {
            for (int i = 0; i < rows; i++)
            {
                if (i != lead)
                {
                    LinearComb(red, i, lead, lead);
                }
            }
        }
        if(red[lead][leadcol] != 0 && red[lead][leadcol] != 1)
        {
            RowMultiply(red, lead, lead);
            for (int i = 0; i < rows; i++)
            {
                if (i != lead)
                {
                    LinearComb(red, i, lead, lead);
                }
            }
        }
        cout<<"\n";
        PrintMatrix(red);
        cout<<"\n";
        lead++;
        leadcol++;
    }
    cout<<"NOW IN SECOND LOOP \n";
    for (int p = 0; p < rows; p++)
    {
        lead = 0;
        leadcol = 0;
        while (lead < rows)
        {
            if (red[lead][leadcol] == 0)
            {
                int counter = 0;
                for (int o = 0; o < red.numcols(); o++)
                {
                    if (red[lead][o] != 0)
                    {
                        counter++;
                    }
                }
                if (counter == 0 && (lead+1 < rows))
                {
                    Swap(red, lead, lead+1);
                }
            }
            lead++;
            leadcol++;
        }
    }
    cout<<"\n";
    PrintMatrix(red);
    cout<<"\n";
    return red;
}
apmatrix<float> Invert(apmatrix<float> in)
{
    int rows = in.numrows();
    int cols = in.numcols();
    apmatrix<float> inv;
    apmatrix<float> temp = augmentIdentity(in);
    PrintMatrix(temp);
    temp = RowReduceMatrix(temp);
    PrintMatrix(temp);
    inv.resize(rows, cols);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < rows; j++)
        {
            inv[i][j] = temp[i][cols+j];
        }
    }
    return inv;
}
void Swap(apmatrix<float> &in, int r1, int r2)
{
    apmatrix<float> temp;
    temp.resize(1, in.numcols());
    for (int i = 0; i < in.numcols(); i++)
    {
        temp[0][i] = in[r1][i];
        in[r1][i] = in[r2][i];
        in[r2][i] = temp[0][i];
        cout<<"Swapping row " <<r1 <<" with row " <<r2 <<"\n";
        PrintMatrix(in);
    }
}
void LinearComb(apmatrix<float> &in, int r, int lead, int leadcol)
{
    float mult = in[r][lead]/in[lead][leadcol];
    for (int i = 0; i < in.numcols(); i++)
    {
        in[r][i] -= mult * in[lead][i];
        cout<<"Subtracting " <<mult <<" times row " <<lead <<" from row " <<r <<"\n";
        PrintMatrix(in);
    }
}
void RowMultiply(apmatrix<float> &in, int lead, int leadcol)
{
    float mult = (1/in[lead][leadcol]);
    for (int i = 0; i < in.numcols(); i++)
    {
        in[lead][i] *= mult;
    }
    cout<<"Multiplying row " <<lead <<" by " <<mult <<"\n";
    PrintMatrix(in);
}
apmatrix<float> createIdentity(int size)
{
    apmatrix<float> iden;
    iden.resize(size,size);
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            iden[i][j] = 0;
        }
        iden[i][i] = 1;
    }
    return iden;
}
apmatrix<float> augmentIdentity(apmatrix<float> in)
{
    int rows = in.numrows();
    int cols = in.numcols();
    if (rows != cols)
    {
        cout<<"unable to augment identity, not passed a square matrix";
        return in;
    }
    else
    {
        apmatrix<float> iden = createIdentity(rows);
        apmatrix<float> augIden;
        augIden.resize(rows,rows+cols);
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                augIden[i][j] = in[i][j];
            }
            for (int k = 0; k < rows; k++)
            {
                augIden[i][cols+k] = iden[i][k];
            }
        }
        return augIden;
    }
}
apmatrix<float> OLS(apmatrix<float> m, apmatrix<float> y)
{
    apmatrix<float> mtmi, mt, mtm, mty, a;
    mt = TransposeMatrix(m);
    //cout<<"MT \n";
    //PrintMatrix(mt);
    mtm = MultiplyMatricies(mt, m);
    //cout<<"MTM \n";
    //PrintMatrix(mtm);
    mty = MultiplyMatricies(mt, y);
    //cout<<"MTY \n";
    //PrintMatrix(mty);
    mtmi = Invert(mtm);
    //cout<<"MTMI \n";
    //PrintMatrix(mtmi);
    a = MultiplyMatricies(mtmi, mty);
    //cout<<"A \n";
    //PrintMatrix(a);
    return a;
}
apmatrix<float> projection(apmatrix<float> thi, apmatrix<float> ontothat)
{
    //cout<<"Start projection method \n";
    apmatrix<float> projection;
    float mult, top, bottom;
    top = dotprodvec(thi, ontothat);
    bottom = pow(vectormagnitude(thi), 2);
    mult = top/bottom;
    projection = ScalarMultiplyMatrix(thi, mult);
    //cout<<"end projection method \n";
    return projection;
}
float dotprodvec(apmatrix<float> u, apmatrix<float> v)
{
    //cout<<"Start dot prod method \n";
    float dotprod = 0;
    for (int i = 0; i < u.numrows(); i++)
    {
        dotprod += u[i][0] * v[i][0];
    }
    //cout<<"end dot prod method \n";
    //cout<<"dotprod = " <<dotprod;
    return dotprod;
}
float vectormagnitude(apmatrix<float> vec)
{
    //cout<<"Start Vector magnitude method \n";
    float temp = 0;
    for (int i = 0; i < vec.numrows(); i++)
    {
        temp += pow(vec[i][0], 2);
    }
    temp = sqrt(temp);
    //cout<<"end vector magnitude method \n";
    //cout<<"vector magnitude = " <<temp;
    return temp;
}
apmatrix<float> ScalarMultiplyMatrix(apmatrix<float> in, float scalar)
{
    //cout<<"Start Scalar Multiply method \n";
    for (int i = 0; i < in.numrows(); i++)
    {
        for (int j = 0; j < in.numcols(); j++)
        {
            in[i][j] = scalar * in[i][j];
        }
    }
    //cout<<"End Scalar Multiply method \n";
    return in;
}
apmatrix<float> VectorSubtract(apmatrix<float> in1, apmatrix<float> in2)
{
    //cout<<"Start Scalar Multiply method \n";
    apmatrix<float> temp;
    temp.resize(in1.numrows(), 1);
    for (int i = 0; i < in1.numrows(); i++)
    {
        temp[i][0] = in1[i][0] - in2[i][0];
    }
    //cout<<"End Scalar Multiply method \n";
    return temp;
}
apmatrix<float> MakeQMatrix(apmatrix<float> a)
{
    int rows = a.numrows();
    int cols = a.numcols();
    apmatrix<float> Q;
    Q.resize(rows, cols);
    apmatrix<float> c[10];
    apmatrix<float> u[10];
    apmatrix<float> umagnorm[10];
    float umag[10];
    //cout<<"Start first loop \n";
    for (int i = 0; i < cols; i++)
    {
        apmatrix<float> temp;
        temp.resize(rows, 1);
        for (int j = 0; j < rows; j++)
        {
            temp[j][0] = a[j][i];
        }
        c[i] = temp;
    }
    //cout<<"end first loop \n";
    //cout<<"Start second loop \n";
    for (int k = 0; k < cols; k++)
    {
        u[k] = c[k];
        if (k != 0)
        {
            for (int p = 0; p < k; p++)
            {
                u[k] = VectorSubtract(u[k], projection(u[p], c[k]));
            }
        }
        umag[k] = vectormagnitude(u[k]);
        umagnorm[k] = ScalarMultiplyMatrix(u[k], 1/umag[k]);
        //PrintMatrix(umagnorm[k]);
    }
    //reconstructing Q from umagnorm vectors
    //cout<<"end second loop \n";
    //cout<<"Start third loop \n";
    for (int t = 0; t < cols; t++)
    {
        for (int y = 0; y < rows; y++)
        {
            float tempentry = (umagnorm[t])[y][0];
            Q[y][t] = tempentry;
        }
    }
    //cout<<"end third loop \n";
    //PrintMatrix(Q);
    return Q;
}
//12    -51    4
//6    167    -68
//-4    24    -41
apmatrix<float> MakeRMatrix(apmatrix<float> A, apmatrix<float> Q)
{
    apmatrix<float> R = MultiplyMatricies(TransposeMatrix(Q), A);
    //PrintMatrix(R);
    return R;
}
apmatrix<float> QR(apmatrix<float> a, int times)
{
    apmatrix<float> A = a;
    apmatrix<float> Q;
    apmatrix<float> R;
    for (int i = 0; i < times; i++)
    {
        cout<<"A \n";
        PrintMatrix(A);
        Q = MakeQMatrix(A);
        cout<<"Q \n";
        PrintMatrix(Q);
        R = MakeRMatrix(A,Q);
        cout<<"R \n";
        PrintMatrix(R);
        A = MultiplyMatricies(R, Q);
    }
    return A;
}
apmatrix<float> geteigenvalues(apmatrix<float> A, int times)
{
    apmatrix<float> eigenvalues;
    A = QR(A, times);
    PrintMatrix(A);
    eigenvalues.resize(A.numrows(), 1);
    for (int i = 0; i < A.numrows(); i++)
    {
        eigenvalues[i][0] = A[i][i];
    }
    return eigenvalues;
}
float pointpointDist(apmatrix<float> A, apmatrix<float> B)
//should take in an n x 1 matrix representing the two points
{
    if (A.numrows() != B.numrows() || A.numcols() != B.numcols())
    {
        cout<<"Error, not matching dimensions" <<"\n";
        return -1;
    }
    else
    {
        float temp = 0;
        for (int i = 0; i < A.numrows(); i++)
        {
            temp += pow((A[i][0] - B[i][0]),2);
        }
        temp = pow(temp,0.5);
        //cout<<"distance: " << temp <<"\n";
        return temp;
    }
}
apmatrix<float> pointlineDist(apmatrix<float> point, apmatrix<float> linePoint, apmatrix<float> lineVec)
{
    //So you create a vector going from point to linepoint called V, which wont be the shortest
    //then you project this onto the direction vector to find the portion of the generated vector V that is unnecessary
    //then you subtract this from the vector V, which results in only the shortest portion needed to span from point to line

    apmatrix<float> v = SubtractMatricies(linePoint, point);
    apmatrix<float> badpart = projection(v, lineVec);
    v = SubtractMatricies(v, badpart);
    
    //MACHINE LEARNING VERSION
    /*
    float init = floatRand();
    float nudge = 0.001;
    float learningRate = 0.05;
    apmatrix<float> testplus;
    apmatrix<float> testminus;
    apmatrix<float> test = AddMatricies(linePoint, ScalarMultiplyMatrix(directionVec, init));
    for(int i = 0; i < 1000000; i++)
    {
        cout<<i <<"\n";
        testplus = AddMatricies(linePoint, ScalarMultiplyMatrix(directionVec, init+nudge));
        testminus = AddMatricies(linePoint, ScalarMultiplyMatrix(directionVec, init-nudge));
        if (pointpointDist(testplus, point) > pointpointDist(test, point) && pointpointDist(testminus, point) > pointpointDist(test, point))
        {
            break;
            
            return test;
        }
        else if(pointpointDist(testplus, point) < pointpointDist(test, point))
        {
            init += learningRate;
            test = AddMatricies(linePoint, ScalarMultiplyMatrix(directionVec, init));
        }
        else if(pointpointDist(testminus, point) < pointpointDist(test, point))
        {
            init -= learningRate;
            test = AddMatricies(linePoint, ScalarMultiplyMatrix(directionVec, init));
        }
    }
    cout<<pointpointDist(test, point) <<"\n";
    return test;
}
float floatRand()
{
    return float(rand()) / (float(RAND_MAX) + 1.0);
}
*/
    return v;
}

