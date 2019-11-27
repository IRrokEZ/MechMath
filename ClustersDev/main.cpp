#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <cstring>
#include <clocale>
#include <cmath>
#include <memory>
#include <sstream>

const double PI = 3.1415926535897932384626433832795;
const double eps = 0.0000001;

using namespace std;

fstream Mylog;

/*

How to write into file:
Mylog << somedata << endl;

*/

int MN(int a, int b){
    return ((a > b) ? b : a);
}

int MX(int a, int b){
    return ((a < b) ? b : a);
}

double GradToRad(double Grad){//convertion of degrees to radians
    return (Grad * PI / 180);
}

void Error(int type){
    switch(type){
        case 1: cout << endl << "Size of vector less than 1" << endl; break;
        case 2: cout << endl << "Size of group less than 1" << endl; break;
        case 3: cout << endl << "Couldn't open file" << endl; break;
        case 4: cout << endl << "Worksize mustn't be less than 1" << endl; break;
        case 5: cout << endl << "Incorrect input. Please try again" << endl; break;
        case 6: cout << endl << "Incorrect value of point' label" << endl; break;
        case 7: cout << endl << "Incorrect file data" << endl; break;
        default: cout << endl << "Unnamed Error" << endl; break;
    }
}

bool InVector(const vector <int> &vect, int elem){//checking of element in vector
    bool rez = false;
    for(unsigned int i = 0; i < vect.size(); i ++){
        if (vect[i] == elem){
            rez = true;
            return true;
        }
    }
    return rez;
}

double sqr(double val){
    return val*val;
}

double Correct(double val){
    if(val < eps){
        return 0.0;
    }
    return val;
}

class Point{//points
    private:
        double x_koord, y_koord;
        int label;
    public:
        Point(){//Default constructor
            this -> x_koord = 0;
            this -> y_koord = 0;
            this -> label = 0;
        }
        Point(double x, double y){// Initialize fields with passed parameters
            this -> x_koord = x;
            this -> y_koord = y;
            this -> label = 0;
        }
        Point(double x, double y, int lb){// Same, but with a label
            this -> x_koord = x;
            this -> y_koord = y;
            this -> label = lb;
        }
        Point(const Point &other){// Copy constructor
            this -> x_koord = other.x_koord;
            this -> y_koord = other.y_koord;
            this -> label = other.label;
        }
        double GetX(){
            return this -> x_koord;
        }
        double GetY(){
            return this -> y_koord;
        }
        int GetLabel(){
            return this -> label;
        }
        void SetX(double x){
            this -> x_koord = x;
        }
        void SetY(double y){
            this -> y_koord = y;
        }
        void SetLabel(int l){
            this -> label = l;
        }
        void AddX(double x){
            this -> x_koord += x;
        }
        void AddY(double y){
            this -> y_koord += y;
        }
        Point& operator = (const Point& other){
            x_koord = other.x_koord;
            y_koord = other.y_koord;
            label = other.label;
            return *this;
        }
};

vector <Point> arrr(const string& filename){//read Point vector from file
    fstream f;
    double x, y;
    int l;
    f.open(filename, ios::in);
    vector <Point> rez;
    while(!f.eof()){
        f >> x >> y >> l;
        rez.push_back(Point(x, y, l));
    }
    return rez;
}

class Control{
    private:
        vector <double> vec1;
        vector <double> vec2;
        vector <Point> group;
        int numpoint1, numpoint2, numpointgroup;
    public:
        Control(){//default constructor
            this -> vec1.resize(0);
            this -> vec2.resize(0);
            this -> group.resize(0);
            this -> numpoint1 = 0;
            this -> numpoint2 = 0;
            this -> numpointgroup = 0;
        }
        void MakeLabel(int lb){
            for(int i = 0; i < this -> numpointgroup; i ++){
                this -> group[i].SetLabel(lb);
            }
        }
        vector <double> CreateNorm(int numberofpoints, double minval, double maxval){//norm gisto
            if(numberofpoints < 1){
                Error(1);
                return vector <double>(); // Empty, default-constructed vector
            }
            vector <double> arr(numberofpoints); // Create a vector of predefined size
            double sum, temp;
            int mmn, mmx, range;
            if(maxval < minval){
                swap(minval, maxval);
            }
            mmn = (int(minval)) * 10;
            mmx = (int(maxval) + 1) * 10;
            range = mmx - mmn + 1;
            for(int i = 0; i < numberofpoints; i ++){
                sum = 0;
                for(int j = 0; j < 1000; j ++){
                    temp = ((rand () % range) + mmn) / 10.0;
                    if((temp < minval) || (temp > maxval)){
                        j --;
                    } else {
                        sum += temp;
                    }
                }
                sum /= 1000;
                if((sum < minval) || (sum > maxval)){
                    i --;
                } else {
                    arr[i] = sum;
                }
            }
            return arr;
        }
        void GenRnd(int numberofpoints, double minval, double maxval){//ravn gisto
            if(numberofpoints < 1){
                Error(1);
            } else {
                vec1.clear();
                this -> numpoint1 = numberofpoints;
                for(int i = 0; i < this -> numpoint1; i ++){
                    this -> vec1.push_back(((rand() % (int(100 * maxval) - int(100 * minval))) + int(100 * minval)) / 100.0);
                }
            }
        }
        void GenNorm(int numberofpoints, double minval, double maxval){//norm generation
            if(numberofpoints < 1){
                Error(1);
            } else {
                this -> vec2.clear();
                this -> numpoint2 = numberofpoints;
                if(maxval < minval){
                    swap(minval, maxval);
                }
                vector <double> arr = this -> CreateNorm(this -> numpoint2, minval, maxval);
                for(int i = 0; i < this -> numpoint2; i ++){
                    vec2.push_back(arr[i]);
                }
            }
        }
        void GenGroup(int numberofpoints, double minvalx, double maxvalx, double minvaly, double maxvaly, int label){//group creating
            if(numberofpoints < 1){
                Error(2);
            } else {
                this -> numpointgroup = numberofpoints;
                this -> group.clear();
                vector <double> arrx = this -> CreateNorm(this -> numpointgroup, minvalx, maxvalx);
                vector <double> arry = this -> CreateNorm(this -> numpointgroup, minvaly, maxvaly);
                for(int i = 0; i < this -> numpointgroup; i ++){
                    group.push_back(Point(arrx[i], arry[i], label));
                }
            }
        }
        void GenGroup(const vector <Point>& values){//group creating
            this -> numpointgroup = values.size();
            this -> group.clear();
            for(int i = 0; i < this -> numpointgroup; i ++){
                this -> group.push_back(Point(values[i]));
            }
        }
        void FileRavn(){//printing ravn vector in file
            fstream text, script;
            script.open("plotravn.plt", ios::out | ios::trunc); //writing script to the file
            script << "width=1" << endl << "bin(x, s) = s*int(x/s) + width/2"
                << endl << "set boxwidth width" << endl << "plot 'ravn.txt' u (bin($1,width)):(1.0) \\"
                << endl << "s f w boxes fs solid 0.5 title 'Ravn Gisto'" << endl;
            script.close();
            text.open("ravn.txt", ios::out | ios::trunc); //writing values to the file
            for(int i = 0; i < this -> numpoint1; i ++){
                text << this -> vec1[i] << endl;
            }
            text.close();
        }
        void FileNorm(){//printing norm vector in file
            fstream text, script;
            script.open("plotnorm.plt", ios::out | ios::trunc);
            script << "width=0.1" << endl << "bin(x, s) = s*int(x/s) + width/2"
                << endl << "set boxwidth width" << endl << "plot 'norm.txt' u (bin($1,width)):(1.0) \\"
                << endl << "s f w boxes fs solid 0.5 title 'Norm Gisto'" << endl;
            script.close();
            text.open("norm.txt", ios::out | ios::trunc);
            for(int i = 0; i < this -> numpoint2; i ++){
                text << this -> vec2[i] << endl;
            }
            text.close();
        }
        void FileGroup(){//printing group in file
            fstream text, script;
            script.open("plotgroup.plt", ios::out | ios::trunc);
            script << "plot 'group.txt'";
            script.close();
            text.open("group.txt", ios::out | ios::trunc); //base file
            for(int i = 0; i < this -> numpointgroup; i ++){
                text << this -> group[i].GetX() << " " << this -> group[i].GetY() << endl;
            }
            text.close();
        }
        vector <double> RetRavn(){//ravn generation cpy
            vector <double> arr(this->numpoint1);
            for(int i = 0; i < this -> numpoint1; i ++){
                arr[i] = this -> vec1[i];
            }
            return arr;
        }
        vector <double> RetNorm(){//norm generation cpy
            vector <double> arr(this -> numpoint2);
            for(int i = 0; i < this -> numpoint2; i ++){
                arr[i] = this -> vec2[i];
            }
            return arr;
        }
        vector <Point> RetGroup(){//group cpy
            vector <Point> arr; // No default ctor for Point, so no resizes
            arr.reserve(this -> numpointgroup); // At least we'll hint about our dataset size
            for(int i = 0; i < this -> numpointgroup; i ++){
                arr.push_back(group[i]);
            }
            return arr;
        }
        void turnNULL(double phi){//rotarion group relative (0;0)
            double newx, newy;
            vector <Point> newvect;
            newvect.reserve(this -> numpointgroup);
            for(int i = 0; i < this -> numpointgroup; i ++){
                newx = this -> group[i].GetX();
                newy = this -> group[i].GetY();
                newvect.push_back(Point(
                    newx * cos(phi) - newy * sin(phi),
                    newx * sin(phi) + newy * cos(phi)));
            }
            this -> group.clear();
            this -> GenGroup(newvect);
        }
        void turnCenter(double phi){//rotarion group relative to group center
            double midx = 0, midy = 0, newx, newy;
            vector <Point> newvect;
            newvect.reserve(this -> numpointgroup);
            for(int i = 0; i < this -> numpointgroup; i ++){
                midx += this -> group[i].GetX();
                midy += this -> group[i].GetY();
            }
            midx /= this -> numpointgroup;
            midy /= this -> numpointgroup;
            for(int i = 0; i < this -> numpointgroup; i ++){
                newx = this -> group[i].GetX() - midx;
                newy = this -> group[i].GetY() - midy;
                newvect.push_back(Point(
                    (newx * cos(phi) - newy * sin(phi)) + midx,
                    (newx * sin(phi) + newy * cos(phi)) + midy));
            }
            this -> group.clear();
            this -> GenGroup(newvect);
        }
        void MoveX(double deltax){//X axis moving
            vector <Point> newvect = this -> RetGroup();
            for(int i = 0; i < this -> numpointgroup; i ++){
                newvect[i].SetX(newvect[i].GetX() + deltax);
            }
            this -> group.clear();
            this -> GenGroup(newvect);
        }
        void MoveY(double deltay){//Y axis moving
            vector <Point> newvect = this -> RetGroup();
            for(int i = 0; i < this -> numpointgroup; i ++){
                newvect[i].SetY(newvect[i].GetY() + deltay);
            }
            this -> group.clear();
            this -> GenGroup(newvect);
        }
        int RetGroupSize(){
            return this -> numpointgroup;
        }
};

class Group{
    private:
        unique_ptr <Control> group;
        int GroupLabel;
    public:
        Group() : group(new Control()),GroupLabel(0) {} // Default constructor, required for vector::resize
        Group(int groupsize, double minx, double miny, double maxx, double maxy, int GL){//constructor
            this -> group.reset(new Control());
            this -> GroupLabel = GL;
            this -> group -> GenGroup(groupsize, minx, maxx, miny, maxy, GL);
        }
        Group(unique_ptr <Control> newgroup){//constructor //it means that label included into Point' vector
            vector <Point> arr = newgroup -> RetGroup();
            this -> GroupLabel = arr[0].GetLabel();
            this -> group.reset(new Control());
            this -> group -> GenGroup(arr);
        }
        Group(unique_ptr <Control> newgroup, int LB){//constructor
            vector <Point> arr = newgroup -> RetGroup();
            this -> GroupLabel = LB;
            this -> group.reset(new Control());
            this -> group -> GenGroup(arr);
            this -> group -> MakeLabel(LB);
        }
        int RetGroupLabel(){
            return this -> GroupLabel;
        }
        unique_ptr <Control> Ret_Field_Group() const {//returning of a group
            unique_ptr <Control> rgroup(new Control());
            vector <Point> arr = this -> group -> RetGroup();
            rgroup -> GenGroup(arr);
            return rgroup;
        }
        Group(const Group &newgroup) {//cpy constructor
            unique_ptr <Control> GR = newgroup.Ret_Field_Group();
            vector <Point> arr = GR -> RetGroup();
            this -> group.reset(new Control());
            this -> GroupLabel = arr[0].GetLabel();
            this -> group -> GenGroup(arr);
        }
        Group(unique_ptr <Group> &newgroup) {//cpy constructor
            unique_ptr <Control> GR = newgroup->Ret_Field_Group();
            vector <Point> arr = GR -> RetGroup();
            this -> group.reset(new Control());
            this -> GroupLabel = arr[0].GetLabel();
            this -> group -> GenGroup(arr);
        }
        void Regrupp(const unique_ptr <Control> &newgroup){ //we modified this object by new Control object
            vector <Point> arr = newgroup -> RetGroup();
            this -> group.reset(new Control());
            this -> group -> GenGroup(arr);
        }
        void Regrupp(const unique_ptr <Control> &newgroup, int LB){ //we modified this object by new Control object
            vector <Point> arr = newgroup -> RetGroup();
            this -> group.reset(new Control());
            this -> group -> GenGroup(arr);
            this -> GroupLabel = LB;
            this -> group -> MakeLabel(LB);
        }
        void RemakeLabel(){
            this -> group -> MakeLabel(this -> GroupLabel);
        }
        void RemakeLabel(int LB){
            this -> GroupLabel = LB;
            this -> group -> MakeLabel(this -> GroupLabel);
        }
};

class Field{
    private:
        vector <Group> field;
        int fieldsize/*number of groups*/, pointer, totalsize/*number of points*/;
        vector <int> arrsz;
        vector <Point> allkoord;
//SV >
        vector <double> allLengths;
        vector <int> label, counters, marks;

        vector <Point> centerKoords, startPoints, save, workPoints, fieldkoord, forClusters, saveData;;

        vector <vector <int>> binary_table;
        vector <vector <double>> RO;

        double bestWeight, maxRO;
        int k, fsize, optimalK, launchNumber/*count each launch this part of program*/;
//SV <
    public:
//SV >
        void BeforeWork(){
            this -> workPoints.resize(0);
            this -> fieldkoord.resize(0);
            this -> allLengths.resize(0); //we will set them after
            this -> centerKoords.resize(0);
            this -> startPoints.resize(0);
            this -> save.resize(0);
            this -> k = 0;
            this -> maxRO = 0;
            this -> optimalK = 0;
            this -> bestWeight = 0;
            this -> fsize = 0;
            this -> label.resize(0);
            this -> binary_table.resize(0);
            this -> RO.resize(0);
            this -> counters.resize(0);
            this -> forClusters.resize(0);
            this -> saveData.resize(0);
        }
//SV <
//SV >
        void BeforeStart(){
            if(saveData.size() > 0){
                this -> workPoints.resize(0);
                for(unsigned int i = 0; i < this -> saveData.size(); i ++){
                    this -> workPoints.push_back(Point(this -> saveData[i]));
                }
                this -> saveData.resize(0);
            }
            for(unsigned int i = 0; i < this -> workPoints.size(); i ++){
                this -> workPoints[i].SetLabel(0);
            }
            this -> fieldkoord = this -> workPoints;
            this -> allLengths.resize(0); //we will set them after
            this -> centerKoords.resize(0);
            this -> startPoints.resize(0);
            this -> save.resize(0);
            this -> k = 0;
            this -> maxRO = 0;
            this -> optimalK = 0;
            this -> bestWeight = 0;
            this -> fsize = this -> workPoints.size();
            this -> label.resize(this -> workPoints.size());
            this -> binary_table.resize(this -> workPoints.size());
            this -> RO.resize(this -> workPoints.size());
            this -> counters.resize(0);
            this -> forClusters.resize(0);
            for(unsigned int i = 0; i < this -> workPoints.size(); i ++){
                this -> binary_table[i].resize(this -> workPoints.size());
                this -> RO[i].resize(this -> workPoints.size());
            }
        }
//SV <
        Field(int numofgroups){//constructor
            this -> BeforeWork();
            if(numofgroups > 0){
                this -> fieldsize = numofgroups;
            } else {
                this -> fieldsize = 0;
            }
            this -> field.resize(0); //Reserve memory for our object
            this -> pointer = 0;
            this -> totalsize = 0;
            this -> arrsz.resize(this -> fieldsize); // int default-constructs to 0
        }
        Field(const vector <Point>& koords){//constructor
            this -> BeforeWork();
            this -> totalsize = koords.size();
            this -> allkoord = koords;
            this -> field.resize(this -> totalsize); //Reserve memory for our object
            this -> arrsz.resize(this -> totalsize); // int default-constructs to 0
            this -> fieldsize = 0;
            this -> pointer = 0;
        }
        Field(const string& filename){//constructor from file
            this -> BeforeWork();
            this -> allkoord = arrr(filename);
            this -> totalsize = allkoord.size();
            this -> field.resize(this -> totalsize); //Reserve memory for our object
            this -> arrsz.resize(this -> totalsize); // int default-constructs to 0
            this -> fieldsize = 0;
            this -> pointer = 0;
        }
        int retTOTSIZE(){
            return this -> totalsize;
        }
        vector <Point> retALLKOORD(){
            return this -> allkoord;
        }
        void AddGroup(const unique_ptr <Group> &newgroup){//adding of a group to the field
            unique_ptr <Control> tc = newgroup -> Ret_Field_Group();
            int gsz = tc -> RetGroupSize();
            this -> arrsz.push_back(gsz);
            unique_ptr <Group> TempGroup(new Group());
            TempGroup -> Regrupp(tc, this -> pointer);
            this -> field.push_back(Group(TempGroup));
            this -> pointer ++;
        }
        void MakeAllKoord(){//list of field coordinates
            int totsize = 0;
            this -> allkoord.clear();
            vector <Point> tmp;
            for(int i = 0; i < this -> fieldsize; i ++){
                tmp.clear();
                unique_ptr <Control> tempcontrol = this -> field[i].Ret_Field_Group();
                tmp = tempcontrol -> RetGroup();
                totsize += tmp.size();
                for(unsigned int j = 0; j < tmp.size(); j ++){
                    this -> allkoord.push_back(Point(tmp[j]));
                }
            }
            this -> totalsize = totsize;
//SV >
            for(unsigned int i = 0; i < this -> allkoord.size(); i ++){
                this -> workPoints.push_back(Point(this -> allkoord[i]));
            }
//SV <
        }
        void ToTxt(){//printing of the field to file
            fstream text;
            string answ = "wait";
            bool addcolor = false;
            while((answ != "yes") && (answ != "no")){
                cout << endl << "Do you want to add color? " << endl << "Please enter 'yes' or 'no'"<<endl;
                cin>>answ;
                if((answ != "yes") && (answ != "no")){
                    cout<<endl<<"Error!"<<endl;
                }
            }
            if(answ == "yes"){
                addcolor = true;
                fstream script;
                script.open("plotField.plt", ios::out | ios::trunc);
                script << "rgb(r,g,b) = 65536 * int(r) + 256 * int(g) + int(b) " << endl;
                script << "plot 'Field.txt' using 1:2:(rgb($3,$4,$5)) with points lc rgb variable";
                script.close();
            }
            int r, g, b;//, ptrgroup = 0, ptr = 0;
            text.open("Field.txt", ios::out | ios::trunc);
            r = rand() % 256;
            g = rand() % 256;
            b = rand() % 256;
            for(int i = 0; i < this -> totalsize; i ++){
                text << this -> allkoord[i].GetX() << " " << this -> allkoord[i].GetY();
                if(addcolor){
                    if((i == 0) || (this -> allkoord[i].GetLabel() != this -> allkoord[i - 1].GetLabel())){
                        r = rand() % 256;
                        g = rand() % 256;
                        b = rand() % 256;
                    }
                    text << " " << r << " " << g << " " << b;
                }
                text << endl;
            }
            text.close();
            text.open("FieldBackup.txt", ios::out | ios::trunc);
            text << "FIELD" << endl << "FIELD" << endl;
            for(int i = 0; i < this -> totalsize; i ++){
                text << this -> allkoord[i].GetX() << " " << this -> allkoord[i].GetY() << " "
                    << this -> allkoord[i].GetLabel();
                if(i != this -> totalsize - 1){
                    text << endl;
                }
            }
            text.close();
        }
        void ToTxt(string answ){//printing of the field to file
            fstream text;
            bool addcolor = false;
            if((answ != "yes") && (answ != "no")){
                Error(5);
            }
            if(answ == "yes"){
                addcolor = true;
                fstream script;
                script.open("plotField.plt", ios::out | ios::trunc);
                script << "rgb(r,g,b) = 65536 * int(r) + 256 * int(g) + int(b) " << endl;
                script << "plot 'Field.txt' using 1:2:(rgb($3,$4,$5)) with points lc rgb variable";
                script.close();
            }
            int r, g, b;//, ptrgroup = 0, ptr = 0;
            text.open("Field.txt", ios::out | ios::trunc);
            r = rand() % 256;
            g = rand() % 256;
            b = rand() % 256;
            for(int i = 0; i < this -> totalsize; i ++){
                text << this -> allkoord[i].GetX() << " " << this -> allkoord[i].GetY();
                if(addcolor){
                    if((i == 0) || (this -> allkoord[i].GetLabel() != this -> allkoord[i - 1].GetLabel())){
                        r = rand() % 256;
                        g = rand() % 256;
                        b = rand() % 256;
                    }
                    text << " " << r << " " << g << " " << b;
                }
                text << endl;
            }
            text.close();
            text.open("FieldBackup.txt", ios::out | ios::trunc);
            text << "FIELD" << endl << "FIELD" << endl;
            for(int i = 0; i < this -> totalsize; i ++){
                text << this -> allkoord[i].GetX() << " " << this -> allkoord[i].GetY() << " "
                    << this -> allkoord[i].GetLabel();
                if(i != this -> totalsize - 1){
                    text << endl;
                }
            }
            text.close();
        }
        void FromTXT(const string& filename){ //we need this function to read new values to field without creating new object
            this -> allkoord = arrr(filename);
            this -> totalsize = allkoord.size();
        }
//SV >
        //WAVE
        void CreateRo(){//creating of matrix of distances
            double x1, y1, x2, y2;
            for(unsigned int i = 0; i < this -> workPoints.size(); i ++){
                this -> label[i] = -1;
                x1 = this -> workPoints[i].GetX();
                y1 = this -> workPoints[i].GetY();
                for(unsigned int j = 0; j < this -> workPoints.size(); j ++){
                    x2 = this -> workPoints[j].GetX();
                    y2 = this -> workPoints[j].GetY();
                    this -> RO[i][j] = sqrt(sqr(x1 - x2) + sqr(y1 - y2));
                    if(this -> RO[i][j] > this -> maxRO){
                        this -> maxRO = this -> RO[i][j];
                    }
                }
            }
        }
        void GenBinary(double threshold){//creating of binary matrix
            for(unsigned int i = 0; i < this -> workPoints.size(); i ++){
                for(unsigned int j = 0; j < this -> workPoints.size(); j ++){
                    if((this -> RO[i][j] > 0) && (this -> RO[i][j] < (threshold + eps))){
                        this -> binary_table[i][j] = 1;
                    } else {
                        this -> binary_table[i][j] = 0;
                    }
                }
            }
        }
        void Wave(){//wave clasters finding
            int tag = -1, controlsumm, add_elem;
            vector <int> positions;
            for(unsigned int q = 0; q < this -> workPoints.size(); q ++){
                positions.clear();
                controlsumm = 0; //If we check all strings, where first elems was included into cluster
                add_elem = 0; //If we add new element to cluster
                for(unsigned int i = q; i < this -> workPoints.size(); i ++){
                    if(this -> label[i] == -1){
                        tag ++;
                        positions.push_back(i);
                        add_elem ++;
                        this -> label[i] = tag;
                        for(unsigned int j = 0; j < this -> workPoints.size(); j ++){
                            if(this -> binary_table[i][j] == 1){
                                positions.push_back(j);
                                add_elem ++;
                            }
                        }
                        i = this -> workPoints.size();//break;
                    }
                }
                while(controlsumm != add_elem){
                    for(unsigned int i = 0; i < positions.size(); i ++){
                        controlsumm ++;
                        for(unsigned int j = 0; j < this -> workPoints.size(); j ++){
                            if(this -> binary_table[positions[i]][j] == 1){
                                if(!(InVector(positions, j))){
                                    add_elem ++;
                                    positions.push_back(j);
                                }
                            }
                        }
                    }
                }
                for(const auto &strnumber : positions){
                    this -> label[strnumber] = tag;
                }
            }
            for(unsigned int i = 0; i < this -> workPoints.size(); i ++){
                this -> workPoints[i].SetLabel(this -> label[i]);
                if(this -> label[i] < 0){
                    Error(6);
                }
            }
        }
        //TXT
        void ToTxtCluster(string name){//printing clasters in file. Files set by user
            vector <vector <int>> rgb;
            int kkk = 0;
            fstream script, myFile, bup;
            string filename = "", scriptname = "plot", backupname;
            for(unsigned int i = 0; i < this -> workPoints.size(); i ++){
                if(this -> workPoints[i].GetLabel() > kkk){
                    kkk = workPoints[i].GetLabel();
                }
            }
            rgb.resize(kkk + 1);
            for(int i = 0; i < (kkk + 1); i ++){
                rgb[i].resize(3);
                rgb[i][0] = rand() % 256;
                rgb[i][1] = rand() % 256;
                rgb[i][2] = rand() % 256;
            }
            for(unsigned i = 0; i < name.size(); i ++){
                if(name[i] == '.'){
                    i = name.size();
                } else {
                    filename += name[i];
                }
            }
            scriptname += filename + ".plt";
            backupname = filename + "Backup.txt";
            filename = name;
            script.open(scriptname, ios::out | ios::trunc);
            script << "rgb(r,g,b) = 65536 * int(r) + 256 * int(g) + int(b) " << endl
                << "plot '" << filename << "' using 1:2:(rgb($3,$4,$5)) with points lc rgb variable";
            script.close();
            myFile.open(filename, ios::out | ios::trunc);
            int l;
            for(unsigned int i = 0; i < this -> workPoints.size(); i ++){
                l = this -> workPoints[i].GetLabel();
                if(l >= 0){
                    myFile << this -> workPoints[i].GetX() << " " << this -> workPoints[i].GetY() << " "
                        << rgb[l][0] << " " << rgb[l][1] << " " << rgb[l][2] << endl;
                }
            }
            myFile.close();
            bup.open(backupname, ios::out | ios::trunc);
            if(name != "Field.txt"){
                bup << "CLUSTER" << endl << "Total clusters: " << kkk + 1 << endl
                    << "Total points: " << this -> workPoints.size() << endl;
                vector <Point> bpoint;
                double midx, midy, leftx, rightx, topy, bottomy;
                int counterPoints;
                bool start;
                for(int i = 0; i <= kkk; i ++){
                    bpoint.resize(0);
                    counterPoints = 0;
                    midx = 0;
                    midy = 0;
                    start = true;
                    for(unsigned int j = 0; j < this -> workPoints.size(); j ++){
                        if(this -> workPoints[j].GetLabel() == i){
                            bpoint.push_back(Point(this -> workPoints[i]));
                            counterPoints ++;
                            midx += this -> workPoints[j].GetX();
                            midy += this -> workPoints[j].GetY();
                            if(start){
                                leftx = this -> workPoints[j].GetX();
                                rightx = this -> workPoints[j].GetX();
                                topy = this -> workPoints[j].GetY();
                                bottomy = this -> workPoints[j].GetY();
                                start = false;
                            } else {
                                if(leftx > this -> workPoints[j].GetX()){
                                    leftx = this -> workPoints[j].GetX();
                                }
                                if(rightx < this -> workPoints[j].GetX()){
                                    rightx = this -> workPoints[j].GetX();
                                }
                                if(topy < this -> workPoints[j].GetY()){
                                    topy = this -> workPoints[j].GetY();
                                }
                                if(bottomy > this -> workPoints[j].GetY()){
                                    bottomy = this -> workPoints[j].GetY();
                                }
                            }
                        }
                    }
                    bup << endl << "id: " << i << endl << "left border: " << leftx << endl << "right border: "
                        << rightx << endl << "top border: " << topy << endl << "bottom border: " << bottomy << endl
                        << "center point: " << midx / counterPoints << "; " << midy / counterPoints << endl
                        << "total points: " << counterPoints << endl << endl;
                }
                bup << "CLUSTER" << endl;
            } else {
                bup << "FIELD" << endl << "FIELD" << endl;
            }
            for(int i = 0; i < this -> fsize; i ++){
                bup << this -> workPoints[i].GetX() << " " << this -> workPoints[i].GetY() << " "
                    << this -> workPoints[i].GetLabel() << endl;
            }
            bup.close();
        }
        void KMeans(int userK){
            int minid, qStart, qEnd, ptr = 0;
            double minVal;
            bool change, esc = false;
            if(userK < 0){
                qStart = 1;
                qEnd = this -> fsize;
            } else {
                qStart = userK;
                qEnd = userK + 1;
            }
            vector <double> MuFunction;
            MuFunction.resize(this -> fsize - 1);
            for(int q = qStart; (q < qEnd && !esc); q ++){
                MuFunction[ptr] = 0;
                this -> k = q;
                vector <vector <double>> RO_selected;
                vector <vector <double>> RO_centers;
                change = true;
                this -> startPoints.clear();
                this -> startPoints.reserve(this -> k);
                RO_selected.resize(this -> k);
                this -> centerKoords.resize(this -> k);
                vector <Point> midPoint;
                vector <int> numberOfPointsInEachCluster(this -> k);
                for(int i = 0; i < this -> k; i ++){ //create first matrix of RO between all points by pairs
                    this -> startPoints.push_back(Point(this -> fieldkoord[i]));
                    midPoint.push_back(Point(0, 0, -5));//-5 for mid points
                    numberOfPointsInEachCluster[i] = 0;
                    RO_selected[i].resize(this -> fsize);
                    for(int j = 0; j < this -> fsize; j++){
                        RO_selected[i][j] = Correct(sqrt(sqr(this -> startPoints[i].GetX() - this -> fieldkoord[j].GetX())
                                                         + sqr(this -> startPoints[i].GetY() - this -> fieldkoord[j].GetY())));
                    }
                }
                for(int j = 0; j < this -> fsize; j ++){//mark each point (which center was the nearest)
                    minVal = RO_selected[0][j];
                    minid = 0;
                    for(int i = 1; i < this -> k; i ++){//find minimal ro from j point to k ros of each cluster
                        if(RO_selected[i][j] < minVal){
                            minVal = RO_selected[i][j];
                            minid = i;
                        }
                    }
                    this -> fieldkoord[j].SetLabel(minid);
                    numberOfPointsInEachCluster[minid] ++;
                }
                for(int i = 0; i < this -> fsize; i ++){ //caclulate mid points
                    int id = this -> fieldkoord[i].GetLabel();
                    midPoint[id].AddX(Correct(this -> fieldkoord[i].GetX() / numberOfPointsInEachCluster[id]));
                    midPoint[id].AddY(Correct(this -> fieldkoord[i].GetY() / numberOfPointsInEachCluster[id]));

                }
                for(int i = 0; i < this -> k; i ++){ //make new koord centers
                    this -> centerKoords[i].SetX(midPoint[i].GetX());
                    this -> centerKoords[i].SetY(midPoint[i].GetY());
                }
                while(change){
                    RO_centers.resize(this -> k);
                    change = false;
                    midPoint.clear();
                    numberOfPointsInEachCluster.clear();
                    midPoint.resize(this -> k);
                    numberOfPointsInEachCluster.resize(this -> k);
                    for(int i = 0; i < this -> k; i ++){ //calc all rors from centers
                        RO_centers[i].resize(this -> fsize);
                        midPoint[i] = Point(0,0,-5);
                        numberOfPointsInEachCluster[i] = 0;
                        for(int j = 0; j < this -> fsize; j ++){ //find new ros from new centers to each point
                            RO_centers[i][j] = Correct(sqrt(sqr(this -> centerKoords[i].GetX() - this -> fieldkoord[j].GetX())
                                                            + sqr(this -> centerKoords[i].GetY() - this -> fieldkoord[j].GetY())));
                        }
                    }
                    for(int j = 0; j < this -> fsize; j ++){
                        minVal = RO_centers[0][j];
                        minid = 0;
                        for(int i = 1; i < this -> k; i ++){
                            if(RO_centers[i][j] < minVal){
                                minVal = RO_centers[i][j];
                                minid = i;
                            }
                        }
                        this -> fieldkoord[j].SetLabel(minid);
                        numberOfPointsInEachCluster[minid] ++;
                    }
                    for(int i = 0; i < this -> fsize; i ++){ //caclulate mid points
                        int id = this -> fieldkoord[i].GetLabel();
                        midPoint[id].AddX(Correct(this -> fieldkoord[i].GetX() / numberOfPointsInEachCluster[id]));
                        midPoint[id].AddY(Correct(this -> fieldkoord[i].GetY() / numberOfPointsInEachCluster[id]));
                    }
                    for(int i = 0; i < this -> k; i ++){//check definitions between new and old center koords
                        double xOld, yOld, xNew, yNew;
                        xOld = this -> centerKoords[i].GetX();
                        yOld = this -> centerKoords[i].GetY();
                        xNew = midPoint[i].GetX();
                        yNew = midPoint[i].GetY();
                        if((abs(xOld - xNew) > eps) || (abs(yOld - yNew) > eps)){
                            change = true;
                            i = this -> k;
                        }
                    }
                    if(change){
                        this -> centerKoords.clear();
                        this -> centerKoords.resize(this -> k);
                        for(int i = 0; i < this -> k; i ++){
                            this -> centerKoords[i].SetX(midPoint[i].GetX());
                            this -> centerKoords[i].SetY(midPoint[i].GetY());
                        }
                    }
                }
                for(int i = 0; i < this -> k; i ++){// here we find optimal k
                    for(int j = 0; j < this -> fsize; j ++){
                        if(this -> fieldkoord[j].GetLabel() == i){
                            for(int p = j + 1; p < this -> fsize; p ++){
                                if(this -> fieldkoord[p].GetLabel() == i){
                                   MuFunction[ptr] += Correct(sqrt(sqr(this -> fieldkoord[j].GetX() - this -> fieldkoord[p].GetX()) +
                                                                    sqr(this -> fieldkoord[j].GetY() - this -> fieldkoord[p].GetY())));
                                }
                            }
                        }
                    }
                }
                for (int i = 0; i < q; i ++){//here we find optimal k
                    for(int j = i + 1; j < q; j ++){
                        MuFunction[ptr] += Correct(sqrt(sqr(this -> centerKoords[i].GetX() - this -> centerKoords[j].GetX()) +
                                                      sqr(this -> centerKoords[i].GetY() - this -> centerKoords[j].GetY())));
                    }
                }
                if((q != 1) && ((qEnd - qStart) > 1)){
                    if(MuFunction[ptr - 1] < MuFunction[ptr]){
                        esc = true;
                        ptr -= 2;
                    }
                }
                ptr ++;
            }
            this -> optimalK = ptr + 1;
            this -> bestWeight = MuFunction[ptr];
            this -> workPoints.clear();
            for(unsigned int i = 0; i < this -> fieldkoord.size(); i ++){
                this -> workPoints.push_back(this -> fieldkoord[i]);
            }
        }
        int RetBestK(){
            return this -> optimalK;
        }
        //SPTR
        void CalculateRo(){//calculate all Ro between all points
            for(unsigned int i = 0; i < this -> workPoints.size(); i ++){
                for(unsigned int j = 0; j < this -> workPoints.size(); j ++){
                    if(i == j){
                        this -> RO[i][j] = -1;
                    } else if(j == 0){
                        this -> RO[i][j] = -1;
                    } else {
                        this -> RO[i][j] = sqrt(sqr(this -> workPoints[i].GetX() - this -> workPoints[j].GetX()) +
                                                sqr(this -> workPoints[i].GetY() - this -> workPoints[j].GetY()));
                    }
                }
            }
        }
        void FillLengths(){
            vector <int> index;
            fstream tree, script, text;
            string filename = "Tree" + to_string(this -> launchNumber) + ".txt", scriptname = "plotTree" + to_string(this -> launchNumber) + ".plt",
                sstr = "TempGisto" + to_string(this -> launchNumber) + ".txt", ssstr = "plotTempGisto" + to_string(this -> launchNumber) + ".plt";
            script.open(scriptname, ios::out | ios::trunc);
            script << "plot '" << filename << "' using 1:2 with lines lc rgb \"black\" lw 2 notitle";
            script.close();
            tree.open(filename, ios::out | ios::trunc);
            index.resize(0);
            index.push_back(0);
            vector <vector <double>> TempRo = this -> RO;
            double minRo, mr;
            int minId, mi, saveIND;
            while(this -> allLengths.size() < this -> workPoints.size() - 1){
                mr = -1;
                for(const auto &id : index){//for each id from index array
                    minRo = TempRo[id][0];
                    for(unsigned int i = 0; i < this -> workPoints.size(); i ++){
                        if((minRo < 0) && (TempRo[id][i] > 0)){
                            minRo = TempRo[id][i];
                            minId = i;
                        } else if ((TempRo[id][i] < minRo) && (TempRo[id][i] > 0)){
                            minRo = TempRo[id][i];
                            minId = i;
                        }
                    }
                    if((mr < 0) || (mr > minRo)){
                        mr = minRo;
                        mi = minId;
                        saveIND = id;
                    }
                }
                tree << this -> workPoints[saveIND].GetX() << " " <<this -> workPoints[saveIND].GetY() << endl
                    << this -> workPoints[mi].GetX() << " " <<this -> workPoints[mi].GetY() << endl << endl;
                this -> allLengths.push_back(mr);
                index.push_back(mi);
                for(unsigned int i = 0; i < this -> workPoints.size(); i ++){
                    TempRo[i][mi] = -1;
                }
            }
            tree.close();
            script.open(ssstr, ios::out | ios::trunc);
            script << "width=0.01" << endl << "bin(x, s) = s*int(x/s) + width/2"
                << endl << "set boxwidth width" << endl << "plot '" << sstr << "' u (bin($1,width)):(1.0) \\"
                << endl << "s f w boxes fs solid 0.5 title 'Porog Gisto'" << endl;
            script.close();
            text.open(sstr, ios::out | ios::trunc); //writing values to the file
            for(unsigned int i = 0; i < (this -> workPoints.size() - 1); i ++){
                text << this -> allLengths[i] << endl;
            }
            text.close();
        }
        //HIERARCHY
        void Hierarchy(int userK){
            double x1, y1, x2, y2, minro, midx, midy;
            vector <Point> newPoints = this -> workPoints;
            int p1, p2/*id of point 1 and point 2*/, numpoints = this -> workPoints.size(), numclusters = userK;
            vector <vector <double>> ro;
            while(numpoints != numclusters){
                ro.resize(newPoints.size());
                for(unsigned int i = 0; i < newPoints.size(); i ++){
                    ro[i].resize(newPoints.size());
                    x1 = this -> workPoints[i].GetX();
                    y1 = this -> workPoints[i].GetY();
                    for(unsigned int j = 0; j < newPoints.size(); j ++){
                        if(i == j){
                            ro[i][j] = -1;
                        } else {
                            x2 = this -> workPoints[j].GetX();
                            y2 = this -> workPoints[j].GetY();
                            ro[i][j] = sqrt(sqr(x1 - x2) + sqr(y1 - y2));
                        }
                    }
                }
                p1 = 0;
                p2 = 0;
                minro = ro[0][0];
                for(unsigned int i = 0; i < ro.size(); i ++){
                    for(unsigned int j = 0; j < ro.size(); j ++){
                        if(((minro < 0) || (minro > ro[i][j])) && (ro[i][j] > 0)){
                            minro = ro[i][j];
                            p1 = i;
                            p2 = j;
                        }
                    }
                }
                x1 = newPoints[p1].GetX();
                x2 = newPoints[p2].GetX();
                y1 = newPoints[p1].GetY();
                y2 = newPoints[p2].GetY();
                midx = (x1 + x2) / 2.0;
                midy = (y1 + y2) / 2.0;
                for(unsigned int i = 0; i < newPoints.size(); i ++){
                    ro[i].erase(ro[i].begin() + MX(p1, p2));
                    ro[i].erase(ro[i].begin() + MN(p1, p2));
                }
                ro.erase(ro.begin() + MX(p1, p2));
                ro.erase(ro.begin() + MN(p1, p2));
                newPoints.erase(newPoints.begin() + MX(p1, p2));
                newPoints.erase(newPoints.begin() + MN(p1, p2));
                newPoints.push_back(Point(midx, midy, 0));
                numpoints --;
            }
            vector <vector <double>>rro;
            rro.resize(newPoints.size());
            for(unsigned int i = 0; i < newPoints.size(); i ++){
                rro[i].resize(this -> workPoints.size());
                x1 = newPoints[i].GetX();
                y1 = newPoints[i].GetY();
                for(unsigned int j = 0; j < this -> workPoints.size(); j ++){
                    x2 = this -> workPoints[j].GetX();
                    y2 = this -> workPoints[j].GetY();
                    rro[i][j] = Correct(sqrt(sqr(x1 - x2) + sqr(y1 - y2)));
                }
            }
            for(unsigned int j = 0; j < this -> workPoints.size(); j ++){
                minro = rro[0][j];
                p1 = 0;
                for(unsigned int i = 0; i < newPoints.size(); i ++){
                    if(rro[i][j] < minro){
                        minro = rro[i][j];
                        p1 = i;
                    }
                }
                this -> workPoints[j].SetLabel(p1);
            }
            this -> k = newPoints.size();
            this -> fsize = this -> workPoints.size();
        }
        //FOREL
        void Forel(double threshold){
            for(unsigned int i = 0; i < this -> workPoints.size(); i ++){
                this -> workPoints[i].SetLabel(-i);
            }
            vector <Point> tempPoints = this -> workPoints;
            int startID, counterPoints, ptr, clustID = 0;
            double midx, midy;
            unique_ptr <Point> center (new Point(tempPoints[0]));
            vector <double> ro;
            vector <int> binary, id;
            bool change;
            while(tempPoints.size() > 0){
                startID = rand () % tempPoints.size();
                ro.resize(tempPoints.size());
                binary.resize(tempPoints.size());
                change = true;
                center.reset(new Point(tempPoints[startID]));
                while(change){
                    midx = 0;
                    midy = 0;
                    counterPoints = 0;
                    for(unsigned int i = 0; i < tempPoints.size(); i ++){
                        ro[i] = sqrt(sqr(center -> GetX() - tempPoints[i].GetX()) + sqr(center -> GetY() - tempPoints[i].GetY()));
                        if(ro[i] < threshold){
                            binary[i] = 1;
                            counterPoints ++;
                            midx += tempPoints[i].GetX();
                            midy += tempPoints[i].GetY();
                        } else {
                            binary[i] = 0;
                        }
                    }
                    midx /= counterPoints;
                    midy /= counterPoints;
                    if((abs(center -> GetX() - midx) < eps) && (abs(center -> GetY() - midy) < eps)){//when center stops (in that iteration)
                        change = false;
                        id.resize(counterPoints);
                        ptr = 0;
                        for(unsigned int i = 0; ((i < tempPoints.size()) && (ptr < counterPoints)); i ++){
                            if(binary[i] == 1){
                                id[ptr] = tempPoints[i].GetLabel();
                                ptr ++;
                            }
                        }
                        ptr = counterPoints - 1;
                        for(unsigned int i = this -> workPoints.size() - 1; /*((i >= 0) && (*/ptr >= 0/*))*/; i --){//making labels
                            if(this -> workPoints[i].GetLabel() == id[ptr]){
                                this -> workPoints[i].SetLabel(clustID);
                                ptr --;
                            }
                        }
                        ptr = counterPoints - 1;
                        for(unsigned int i = tempPoints.size() - 1; /*((i >= 0) && (*/ptr >= 0/*))*/; i --){//delete found cluster
                            if(tempPoints[i].GetLabel() == id[ptr]){
                                tempPoints.erase(tempPoints.begin() + i);
                                ptr --;
                            }
                        }
                        clustID ++;
                    } else {
                        center.reset(new Point(midx, midy, -5));
                    }
                }
            }
        }
        //RECOVERY
        void GetFromFile(){
            string filename, specialCommand, tempstr;
            cout << endl << "Please enter backup filename" << endl;
            cin >> filename;
            fstream bup;
            bool read;
            double x, y;
            int l;
            vector <Point> tempvec;
            tempvec.resize(0);
            bup.open(filename, ios::in);
            if(!bup.is_open()){
                cout << endl << "File open error" << endl;
            } else {
                read = true;
                getline(bup, specialCommand);
                if((specialCommand != "CLUSTER") && (specialCommand != "FIELD")){
                    Error(7);
                } else {
                    while(!bup.eof()){
                        if(read){
                            getline(bup, tempstr);
                            if(tempstr == specialCommand){
                                read = false;
                            } else {
                                cout << tempstr << endl;
                            }
                        } else {
                            bup >> x >> y >> l;
                            tempvec.push_back(Point(x, y, l));
                        }
                    }
                }
            }
            bup.close();
            this -> workPoints = tempvec;
            this -> BeforeStart();
        }
        //DBSCAN
        void CalcRo(){
            double x1, y1, x2, y2;
            this -> RO.resize(this -> workPoints.size());
            for(unsigned int i = 0; i < this -> workPoints.size(); i ++){
                this -> RO[i].resize(this -> workPoints.size());
                x1 = this -> workPoints[i].GetX();
                y1 = this -> workPoints[i].GetY();
                for(unsigned int j = 0; j < this -> workPoints.size(); j ++){
                    x2 = this -> workPoints[j].GetX();
                    y2 = this -> workPoints[j].GetY();
                    this -> RO[i][j] = sqrt(sqr(x1 - x2) + sqr(y1 - y2));
                }
            }
        }
        vector <int> CountNearPoints(vector <vector <double>> RoMatrix, double radius){
            vector <int> ctr;
            ctr.resize(this -> workPoints.size());
            for(unsigned int i = 0; i < this -> workPoints.size(); i ++){
                ctr[i] = 0;
                for(unsigned int j = 0; j < this -> workPoints.size(); j ++){
                    if(RoMatrix[i][j] < (radius + eps)){
                        ctr[i] ++;
                    }
                }
            }
            return ctr;
        }
        vector <int> CreateMarks(vector <vector <double>> ro, double radius, vector <int> cnt, int EntNum){
            vector <int> marker;
            marker.resize(this -> workPoints.size());
            for(unsigned int i = 0; i < this -> workPoints.size(); i ++){
                if(cnt[i] > EntNum){
                    marker[i] = 1;
                } else {
                    marker[i] = -1;
                }
            }
            for(unsigned int i = 0; i < this -> workPoints.size(); i ++){
                if(marker[i] == -1){
                    for(unsigned int j = 0; j < this -> workPoints.size(); j ++){
                        if((ro[i][j] < (radius + eps)) && (marker[j] == 1)){
                            marker[i] = 0;
                        }
                    }
                }
            }
            return marker;
        }
        void DelPoints(double radius, int EntryNum){
            this -> fieldkoord.resize(0);
            for(unsigned int i = 0; i < this -> workPoints.size(); i ++){
                this -> fieldkoord.push_back(Point(this -> workPoints[i]));
            }
            this -> CalcRo();
            this -> counters = CountNearPoints(this -> RO, radius);
            this -> marks = CreateMarks(this -> RO, radius, this -> counters, EntryNum);
            for(int i = this -> fieldkoord.size() - 1; i >= 0; i --){
                if(this -> marks[i] == -1){
                    this -> fieldkoord.erase(this -> fieldkoord.begin() + i);
                    this -> counters.erase(this -> counters.begin() + i);
                    for(unsigned int j = 0; j < this -> fieldkoord.size(); j ++){
                        this -> RO[j].erase(this -> RO[j].begin() + i);
                    }
                    this -> RO.erase(this -> RO.begin() + i);
                    this -> marks.erase(this -> marks.begin() + i);
                    this -> workPoints[i].SetLabel(-9999);
                }
            }
        }
        void DBSCANAlghorithm(double radius, int numpointsincircle){
            this -> BeforeStart();
            saveData.resize(0);
            for(unsigned int i = 0; i < this -> workPoints.size(); i ++){
                saveData.push_back(Point(this -> workPoints[i]));
            }
            this -> DelPoints(radius, numpointsincircle);
            this -> workPoints.resize(0);
            for(unsigned int i = 0; i < this -> fieldkoord.size(); i ++){
                this -> workPoints.push_back(Point(this -> fieldkoord[i]));
            }
            this -> BeforeStart();
            this -> CreateRo();
            this -> GenBinary(radius + eps);
            this -> Wave();
        }
//SV <
};

class Cluster{
    private:
        unique_ptr <Field> MyField;
        int launchNumber;
    public:
        Cluster(){//default constructor
            fstream launcher;
            launcher.open("launch.log",  ios::out | ios::trunc);
            launcher << 1;
            launcher.close();
            this -> MyField.reset(new Field(0));
        }
        Cluster(Field *mf){//constructor
            fstream launcher;
            launcher.open("launch.log", ios::in | ios::out);
            launcher >> this -> launchNumber;
            launcher.close();

            this -> MyField.reset(new Field(mf -> retALLKOORD()));
        }
        //RUN WAVE
        void FindByWaveAlgorithm(){
            this -> MyField -> BeforeStart();
            this -> MyField -> CreateRo();
            double threshold = -1;
            string filename = "Wave" + to_string(this -> launchNumber) + ".txt", accept = "zzz";
            fstream launcher;
            this -> launchNumber ++;
            launcher.open("launch.log", ios::in | ios::out | ios::trunc);
            launcher << this -> launchNumber;
            launcher.close();
            while(threshold < 0){
                cout << endl << "Please enter your threshold value" << endl;
                cin >> threshold;
                if(threshold < 0){
                    Error(5);
                }
            }
            this -> MyField -> GenBinary(threshold);
            this -> MyField -> Wave();
            accept = "zzz";
            while((accept != "yes") && (accept != "no")){
                cout << endl << "Do you want to set your file name? ('yes'/'no')" << endl;
                cin >> accept;
                if((accept != "yes") && (accept != "no")){
                    Error(5);
                }
            }
            if(accept == "yes"){
                cout << endl << "Please enter your filename" << endl;
                cin >> filename;
            }
            this -> MyField -> ToTxtCluster(filename);
        }
        //RUN KMEANS
        void FindByKMeansAlgorithm(){
            this -> MyField -> BeforeStart();
            int userNumberOfClusters = -100;
            while((userNumberOfClusters != -1) && (userNumberOfClusters < 0)){
                cout << endl << "Please enter your number of clusters or '-1' to find optimal number" << endl;
                cin >> userNumberOfClusters;
                if(userNumberOfClusters > int(this -> MyField -> retTOTSIZE())){
                    userNumberOfClusters = -1000;
                    cout << endl << "Too big number" << endl;
                }
                if((userNumberOfClusters != -1) && (userNumberOfClusters < 0)){
                    Error(5);
                }
            }
            this -> MyField -> KMeans(userNumberOfClusters);
            if(userNumberOfClusters == -1){
                cout << endl << "Optimal number of clusters is " << this -> MyField -> RetBestK() << endl
                << "You could rerun k-means algorithm with that value" << endl;
            } else {
                string filename = "K-Means" + to_string(this -> launchNumber) + ".txt", accept = "zzz";
                fstream launcher;
                launcher.open("launch.log", ios::in | ios::out | ios::trunc);
                this -> launchNumber ++;
                launcher << this -> launchNumber;
                launcher.close();
                while((accept != "yes") && (accept != "no")){
                    cout << endl << "Do you want to set your file name? ('yes'/'no')" << endl;
                    cin >> accept;
                    if((accept != "yes") && (accept != "no")){
                        Error(5);
                    }
                }
                if(accept == "yes"){
                    cout << endl << "Please enter your filename" << endl;
                    cin >> filename;
                }
                this -> MyField -> ToTxtCluster(filename);
            }
        }
        //RUN SPTR
        void RunSpainningTreeAlgorithm(){
            string scn1 = "TempGisto" + to_string(this -> launchNumber);
            fstream launch;
            this -> MyField -> BeforeStart();
            this -> MyField -> CalculateRo();
            this -> MyField -> FillLengths();
            cout << endl << "You could plot '" << scn1 << "' and choose correct threshold value, then rerun wave algorithm" << endl;
            this -> launchNumber ++;
            launch.open("launch.log", ios::in | ios::out | ios::trunc);
            launch << this -> launchNumber;
            launch.close();
            this -> FindByWaveAlgorithm();
        }
        //RUN HIERARCHY
        void FindByHierarchyAlgorithm(){
            this -> MyField -> BeforeStart();
            int userNumberOfClusters = -100;
            while((userNumberOfClusters != -1) && (userNumberOfClusters < 0)){
                cout << endl << "Please enter your number of clusters" << endl;
                cin >> userNumberOfClusters;
                if(userNumberOfClusters > int(this -> MyField -> retTOTSIZE())){
                    userNumberOfClusters = -1000;
                }
                if(userNumberOfClusters < 0){
                    Error(5);
                }
            }
            this -> MyField -> Hierarchy(userNumberOfClusters);
            string filename = "Hierarchy" + to_string(this -> launchNumber) + ".txt", accept = "zzz";
            fstream launcher;
            launcher.open("launch.log", ios::in | ios::out | ios::trunc);
            this -> launchNumber ++;
            launcher << this -> launchNumber;
            launcher.close();
            while((accept != "yes") && (accept != "no")){
                cout << endl << "Do you want to set your file name? ('yes'/'no')" << endl;
                cin >> accept;
                if((accept != "yes") && (accept != "no")){
                    Error(5);
                }
            }
            if(accept == "yes"){
                cout << endl << "Please enter your filename" << endl;
                cin >> filename;
            }
            this -> MyField -> ToTxtCluster(filename);
        }
        //RUN FOREL
        void FindByForelAlghorithm(){
            this -> MyField -> BeforeStart();
            double threshold = -1;
            while(threshold < 0){
                cout << endl << "Please enter your threshold" << endl;
                cin >> threshold;
                if(threshold < 0){
                    Error(5);
                }
            }
            this -> MyField -> Forel(threshold);
            string filename = "Forel" + to_string(this -> launchNumber) + ".txt", accept = "zzz";
            fstream launcher;
            this -> launchNumber ++;
            launcher.open("launch.log", ios::in | ios::out | ios::trunc);
            launcher << this -> launchNumber;
            launcher.close();
            while((accept != "yes") && (accept != "no")){
                cout << endl << "Do you want to set your file name? ('yes'/'no')" << endl;
                cin >> accept;
                if((accept != "yes") && (accept != "no")){
                    Error(5);
                }
            }
            if(accept == "yes"){
                cout << endl << "Please enter your filename" << endl;
                cin >> filename;
            }
            this -> MyField -> ToTxtCluster(filename);
        }
        //RUN DBSCAN
        void FindByDBSCANAlghorithm(){
            double radius = -1;
            while(radius < 0){
                cout << endl << "Please enter your radius" << endl;
                cin >> radius;
                if(radius < 0){
                    Error(5);
                }
            }
            int NumNearPoints = -1;
            while(NumNearPoints < 0){
                cout << endl << "Please enter number of points in circle" << endl;
                cin >> NumNearPoints;
                if(NumNearPoints < 0){
                    Error(5);
                }
            }
            this -> MyField -> DBSCANAlghorithm(radius, NumNearPoints);
            string filename = "DBSCAN" + to_string(this -> launchNumber) + ".txt", accept = "zzz";
            fstream launcher;
            this -> launchNumber ++;
            launcher.open("launch.log", ios::in | ios::out | ios::trunc);
            launcher << this -> launchNumber;
            launcher.close();
            while((accept != "yes") && (accept != "no")){
                cout << endl << "Do you want to set your file name? ('yes'/'no')" << endl;
                cin >> accept;
                if((accept != "yes") && (accept != "no")){
                    Error(5);
                }
            }
            if(accept == "yes"){
                cout << endl << "Please enter your filename" << endl;
                cin >> filename;
            }
            this -> MyField -> ToTxtCluster(filename);
        }
};

class Interface{
    private:
        int sz;
    public:
        string command;

        Interface(){//constructor
            this -> sz = 1000;
        }
        void PrintHelp(){//help printing in console
            cout << endl << " Help:" << endl;
            cout << " 'genrnd' for DEMO ravn generation of the vector" << endl;
            cout << " 'gennorm' for DEMO norm generation of the vector" << endl;
            cout << " 'genfield' to make field of groups" << endl;
            cout << " 'gengroup' to make one DEMO group" << endl;
            cout << " 'moveX' to move X DEMO" << endl;
            cout << " 'moveY' to move Y DEMO" << endl;
            cout << " 'rollN' to turn group(0;0) DEMO" << endl;
            cout << " 'rollC' to turn group (center) DEMO" << endl;
            cout << " 'exit' to exit" << endl << endl;
            cout << " 'setsize' to set new basic size" << endl << endl;
        }
         void run(){//method of using of the commands by the user
            string command = "start";
            int trg;
            PrintHelp();
            while(command != "exit"){
                trg = -1;
                cout << endl << " Please enter the command" << endl;
                cin >> command;
                if((command == "setsize")){//generation of ravn vector (DEMO)
                    int newsz = 0;
                    while(newsz < 1){
                        cout << endl << "Please enter size (number of points)" << endl;
                        cin >> newsz;
                        if(newsz < 1){
                            Error(4);
                        }
                    }
                    this -> sz = newsz;
                    trg = 0;
                }
                if((command == "genrnd")){//generation of ravn vector (DEMO)
                    double mn, mx;
                    cout << endl << " Please enter min and max value of random" << endl;
                    cin >> mn >> mx;
                    unique_ptr <Control> RND(new Control());
                    RND -> GenRnd(this -> sz, mn, mx);
                    RND -> FileRavn();
                    trg = 0;
                }
                else if(command == "gennorm"){//generation of norm vector (DEMO)
                    double mn, mx;
                    cout << endl << " Please enter min and max value of random" << endl;
                    cin >> mn >> mx;
                    unique_ptr <Control> NORM(new Control());
                    NORM -> GenNorm(this -> sz, mn, mx);
                    NORM -> FileNorm();
                    trg = 0;
                }
                else if(command == "gengroup"){//generation of a group of points (DEMO)
                    double mnx, mxx, mny, mxy;
                    cout << endl << " Please enter min and max value of random (x)" << endl;
                    cin >> mnx >> mxx;
                    cout << endl << " Please enter min and max value of random (y)" << endl;
                    cin >> mny >> mxy;
                    unique_ptr <Control> GROUP(new Control());
                    GROUP -> GenGroup(this -> sz, mnx, mxx, mny, mxy, 0);
                    GROUP -> FileGroup();
                    trg = 0;
                }
                else if(command == "rollN"){//turning of the group (0;0) (DEMO)
                    double phi, mnx, mny, mxx, mxy;
                    vector <Point> gr;
                    cout << endl << "Please create new group to show this function" << endl << endl;
                    cout << endl << " Please enter min and max value of random (x)" << endl;
                    cin >> mnx >> mxx;
                    cout << endl << " Please enter min and max value of random (y)" << endl;
                    cin >> mny >> mxy;
                    unique_ptr <Control> work (new Control());
                    work -> GenGroup(this -> sz, mnx, mxx, mny, mxy, 0); //base group
                    fstream rotated;
                    rotated.open("rotateN.txt", ios::out | ios::trunc);
                    if(!rotated.is_open()){
                        Error(3);
                    }
                    work -> FileGroup();
                    cout << endl << " Please enter angle" << endl;
                    cin >> phi;
                    phi = GradToRad(phi);
                    work -> turnNULL(phi);
                    gr = work -> RetGroup();
                    for(int i = 0; i < this -> sz; i ++){
                        rotated << gr[i].GetX() << " ";
                        rotated << gr[i].GetY() << "\n";
                    }
                    rotated.close();
                    trg = 0;
                }
                else if(command == "rollC"){//turning of the group (Center) (DEMO)
                    double phi, mnx, mny, mxx, mxy;
                    vector <Point> gr;
                    cout << endl << "Please create new group to show this function" << endl;
                    cout << endl << " Please enter min and max value of random (x)" << endl;
                    cin >> mnx >> mxx;
                    cout << endl << " Please enter min and max value of random (y)" << endl;
                    cin >> mny >> mxy;
                    unique_ptr <Control> work (new Control());
                    work -> GenGroup(this -> sz, mnx, mxx, mny, mxy, 0); //base group
                    fstream rotated;
                    rotated.open("rotS.txt", ios::out | ios::trunc);
                    if(!rotated.is_open()){
                        Error(3);
                    }
                    work -> FileGroup();
                    cout << endl << " Please enter angle" << endl;
                    cin >> phi;
                    phi = GradToRad(phi);
                    work -> turnCenter(phi);
                    gr = work -> RetGroup();
                    for(int i = 0; i < this -> sz; i ++){
                        rotated << gr[i].GetX() << " ";
                        rotated << gr[i].GetY() << "\n";
                    }
                    rotated.close();
                    trg = 0;
                }
                else if(command == "moveX"){//moving X (DEMO)
                    double delx, mnx, mny, mxx, mxy;
                    vector <Point> gr;
                    cout << endl << " Please create new group to show this function" << endl;
                    cout << endl << " Please enter min and max value of random (x)" << endl;
                    cin >> mnx >> mxx;
                    cout << endl << " Please enter min and max value of random (y)" << endl;
                    cin >> mny >> mxy;
                    unique_ptr <Control> work (new Control());
                    work -> GenGroup(this -> sz, mnx, mxx, mny, mxy, 0); //base group
                    fstream moved;
                    moved.open("moveX.txt", ios::out | ios::trunc);
                    if(!moved.is_open()){
                        Error(3);
                    }
                    work -> FileGroup();
                    cout << endl << " Please enter delta X" << endl;
                    cin >> delx;
                    work -> MoveX(delx);
                    gr = work -> RetGroup();
                    for(int i = 0; i < this -> sz; i ++){
                        moved << gr[i].GetX() << " ";
                        moved << gr[i].GetY() << "\n";
                    }
                    moved.close();
                    trg = 0;
                }
                else if(command == "moveY"){//moving Y (DEMO)
                    double dely, mnx, mny, mxx, mxy;
                    vector <Point> gr;
                    cout << endl << " Please create new group to show this function" << endl;
                    cout << endl << " Please enter min and max value of random (x)" << endl;
                    cin >> mnx >> mxx;
                    cout << endl << " Please enter min and max value of random (y)" << endl;
                    cin >> mny >> mxy;
                    unique_ptr <Control> work (new Control());
                    work -> GenGroup(this -> sz, mnx, mxx, mny, mxy, 0); //base group
                    fstream moved;
                    moved.open("moveY.txt", ios::out | ios::trunc);
                    if(!moved.is_open()){
                        Error(3);
                    }
                    work -> FileGroup();
                    cout << endl << " Please enter delta Y" << endl;
                    cin >> dely;
                    work -> MoveY(dely);
                    gr = work -> RetGroup();
                    for(int i = 0; i < this -> sz; i ++){
                        moved << gr[i].GetX() << " ";
                        moved << gr[i].GetY() << "\n";
                    }
                    moved.close();
                    trg = 0;
                }
                else if(command == "genfield"){//field creating
                    int numofgroups = 0, tr;
                    double mnx, mny, mxx, mxy;
                    string com = "rrr";
                    unique_ptr <Field> MyField;
                    unique_ptr <Cluster> MyCluster(new Cluster());
                    unique_ptr <Group> temp;
                    while((com != "get") && (com!="work")){
                        cout << endl << "Please enter 'get' to read field from file" << endl;
                        cout << endl << "Please enter 'work' to work in basik mode" << endl;
                        cin >> com;
                        if((com != "get") && (com!="work")){
                            Error(5);//incorrect value
                        }
                    }
                    if(com == "work"){
                        while(numofgroups < 1){
                            cout << endl << " Please enter number of groups in field (>=1)" << endl;
                            cin >> numofgroups;
                            if(numofgroups < 1){
                                cout << endl << " Error! Please enter correctly value" << endl;
                            }
                        }
                        MyField.reset(new Field(numofgroups));
                        for(int i = 0; i < numofgroups; i ++){
                            com = "start";
                            cout << endl << " Please enter min and max value of random (x)" << endl;
                            cin >> mnx >> mxx;
                            cout << endl << " Please enter min and max value of random (y)" << endl;
                            cin >> mny >> mxy;
                            temp.reset(new Group(this -> sz, mnx, mny, mxx, mxy, i));
                            while(com != "esc"){
                                tr = -1;
                                cout << endl << " Enter 'RN' to turn group(0;0)" << endl;
                                cout << endl << " Enter 'RC' to turn group (center)" << endl;
                                cout << endl << " Enter 'MX' to move X" << endl;
                                cout << endl << " Enter 'MY' to move Y" << endl;
                                cout << endl << " Enter 'esc' to finish field creating" << endl;
                                cin >> com;
                                if(com == "RN"){//turning of the group (0;0)
                                    double alf;
                                    cout << endl << " Please enter angle" << endl;
                                    cin >> alf;
                                    alf = GradToRad(alf);
                                    unique_ptr <Control> work = temp -> Ret_Field_Group();
                                    work -> turnNULL(alf);
                                    temp -> Regrupp(work);
                                    tr = 0;
                                }
                                else if(com == "RC"){//turning of the group (Center)
                                    double alf;
                                    cout << endl << " Please enter angle" << endl;
                                    cin >> alf;
                                    alf = GradToRad(alf);
                                    unique_ptr <Control> work = temp -> Ret_Field_Group();
                                    work -> turnCenter(alf);
                                    temp -> Regrupp(work);
                                    tr = 0;
                                }
                                else if(com == "MX"){//moving X
                                    double dx;
                                    unique_ptr <Control> work = temp -> Ret_Field_Group();
                                    cout << endl << " Enter delta x" << endl;
                                    cin >> dx;
                                    work -> MoveX(dx);
                                    temp -> Regrupp(work);
                                    tr = 0;
                                }
                                else if(com == "MY"){//moving Y
                                    double dy;
                                    unique_ptr <Control> work = temp -> Ret_Field_Group();
                                    cout << endl << " Enter delta y" << endl;
                                    cin >> dy;
                                    work -> MoveY(dy);
                                    temp -> Regrupp(work);
                                    tr = 0;
                                }
                                if(tr == -1){
                                    if(com != "esc"){
                                        Error(5); //incorrect value
                                        com = "start";
                                    }
                                }
                            }
                            MyField -> AddGroup(temp);
                        }
                    }
                    else if (com == "get"){
                        MyCluster -> GetFromFile();
                        MyCluster -> ToTxtCluster("Field.txt");
                    }
                    if(com != "get"){
                        MyField -> MakeAllKoord();
                        MyField -> ToTxt();
                        MyCluster.reset(new Cluster(MyField.get()));
                    }
                    string q = "NO";
                    while((q != "yes") && (q != "no")){
                        cout << endl << "Do you want to find clusters? Enter 'yes' or 'no'" << endl;
                        cin >> q;
                        if((q != "yes") && (q != "no")){
                            Error(5);
                        }
                    }
                    if(q == "no"){
                        q = "esc";
                    }
                    while(q != "esc"){
                        cout << endl << "To Find by Wave enter 'wave'" << endl << "To Find by K-means enter 'km'" << endl
                        << "To Spainning Tree enter 'sptr'" << endl << "To Hierarchy enter 'ie'" << endl
                        << "To Forel alghorithm enter 'fish'" << endl << "To DBSCAN enter 'dbs'" << endl
                        << "To Exit enter 'esc'" << endl;
                        cin >> q;
                        if(q == "wave"){
                            MyCluster -> FindByWaveAlgorithm();
                        } else if (q == "km"){
                            MyCluster -> FindByKMeansAlgorithm();
                        } else if (q == "sptr"){
                            MyCluster -> RunSpainningTreeAlgorithm();
                        } else if (q == "ie"){
                            MyCluster -> FindByHierarchyAlgorithm();
                        } else if (q == "fish"){
                            MyCluster -> FindByForelAlghorithm();
                        } else if (q == "dbs"){
                            MyCluster -> FindByDBSCANAlghorithm();
                        } else if (q != "esc"){
                            Error(5);
                        }
                    }
                    trg = 0;
                }
                if(trg == -1){
                    if(command == "exit"){
                        break;
                    } else{
                        cout << " Error! Enter correct command" << endl;
                    }
                }
                else if(trg == 0){
                    cout << endl << " Please enter the command" << endl;
                    PrintHelp();
                    getline(cin, command);
                    trg = -1;
                }
            }
        }
};

class InterfaceSTR{
    private:
        int sz;
        vector <string> comlist;
    public:
        string command;

        InterfaceSTR(){//constructor
            this -> sz = 1000;
            string fname, tempstr;
            this -> comlist.resize(0);
            cout << endl << "Please enter master filename" << endl;
            getline(cin, fname);
            //cin >> fname;
            fstream master;
            master.open(fname, ios::in);
            if(!master.is_open()){
                cout << endl << "File open error" << endl;
            } else {
                while(!master.eof()){
                    getline(master, fname);
                    //cout << endl << fname;
                    if(fname != ""){
                        comlist.push_back(string(fname));
                    }
                }
            }
        }
        void run(){//method of using of the commands by the user
            string command = "start";
            stringstream strm;
            int trg, pos = 0;
            while(command != "exit"){
                trg = -1;
                command = this -> comlist[pos];
                pos ++;
                if((command == "setsize")){//generation of ravn vector (DEMO)
                    int newsz = 0;
                    while(newsz < 1){
                        strm << this -> comlist[pos];
                        pos ++;
                        strm >> newsz;
                        strm.clear();
                        if(newsz < 1){
                            Error(4);
                        }
                    }
                    this -> sz = newsz;
                    trg = 0;
                }
                if((command == "genrnd")){//generation of ravn vector (DEMO)
                    double mn, mx;
                    strm << this -> comlist[pos]; //we hope that file has structure <mn> <mx> in 1 line
                    pos ++;
                    strm >> mn >> mx;
                    strm.clear();
                    unique_ptr <Control> RND(new Control());
                    RND -> GenRnd(this -> sz, mn, mx);
                    RND -> FileRavn();
                    trg = 0;
                }
                else if(command == "gennorm"){//generation of norm vector (DEMO)
                    double mn, mx;
                    strm << this -> comlist[pos]; //we hope that file has structure <mn> <mx> in 1 line
                    pos ++;
                    strm >> mn >> mx;
                    strm.clear();
                    unique_ptr <Control> NORM(new Control());
                    NORM -> GenNorm(this -> sz, mn, mx);
                    NORM -> FileNorm();
                    trg = 0;
                }
                else if(command == "gengroup"){//generation of a group of points (DEMO)
                    double mnx, mxx, mny, mxy;
                    strm << this -> comlist[pos]; //we hope that file has structure <mn> <mx> in 1 line
                    pos ++;
                    strm >> mnx >> mxx;
                    strm.clear();
                    strm << this -> comlist[pos]; //we hope that file has structure <mn> <mx> in 1 line
                    pos ++;
                    strm >> mny >> mxy;
                    strm.clear();
                    unique_ptr <Control> GROUP(new Control());
                    GROUP -> GenGroup(this -> sz, mnx, mxx, mny, mxy, 0);
                    GROUP -> FileGroup();
                    trg = 0;
                }
                else if(command == "rollN"){//turning of the group (0;0) (DEMO)
                    double phi, mnx, mny, mxx, mxy;
                    vector <Point> gr;
                    strm << this -> comlist[pos]; //we hope that file has structure <mn> <mx> in 1 line
                    pos ++;
                    strm >> mnx >> mxx;
                    strm.clear();
                    strm << this -> comlist[pos]; //we hope that file has structure <mn> <mx> in 1 line
                    pos ++;
                    strm >> mny >> mxy;
                    strm.clear();
                    unique_ptr <Control> work (new Control());
                    work -> GenGroup(this -> sz, mnx, mxx, mny, mxy, 0); //base group
                    fstream rotated;
                    rotated.open("rotateN.txt", ios::out | ios::trunc);
                    if(!rotated.is_open()){
                        Error(3);
                    }
                    work -> FileGroup();
                    strm << this -> comlist[pos];
                    pos ++;
                    strm >> phi;
                    strm.clear();
                    phi = GradToRad(phi);
                    work -> turnNULL(phi);
                    gr = work -> RetGroup();
                    for(int i = 0; i < this -> sz; i ++){
                        rotated << gr[i].GetX() << " ";
                        rotated << gr[i].GetY() << "\n";
                    }
                    rotated.close();
                    trg = 0;
                }
                else if(command == "rollC"){//turning of the group (Center) (DEMO)
                    double phi, mnx, mny, mxx, mxy;
                    vector <Point> gr;
                    strm << this -> comlist[pos]; //we hope that file has structure <mn> <mx> in 1 line
                    pos ++;
                    strm >> mnx >> mxx;
                    strm.clear();
                    strm << this -> comlist[pos]; //we hope that file has structure <mn> <mx> in 1 line
                    pos ++;
                    strm >> mny >> mxy;
                    strm.clear();
                    unique_ptr <Control> work (new Control());
                    work -> GenGroup(this -> sz, mnx, mxx, mny, mxy, 0); //base group
                    fstream rotated;
                    rotated.open("rotS.txt", ios::out | ios::trunc);
                    if(!rotated.is_open()){
                        Error(3);
                    }
                    work -> FileGroup();
                    strm << this -> comlist[pos];
                    pos ++;
                    strm >> phi;
                    strm.clear();
                    phi = GradToRad(phi);
                    work -> turnCenter(phi);
                    gr = work -> RetGroup();
                    for(int i = 0; i < this -> sz; i ++){
                        rotated << gr[i].GetX() << " ";
                        rotated << gr[i].GetY() << "\n";
                    }
                    rotated.close();
                    trg = 0;
                }
                else if(command == "moveX"){//moving X (DEMO)
                    double delx, mnx, mny, mxx, mxy;
                    vector <Point> gr;
                    strm << this -> comlist[pos]; //we hope that file has structure <mn> <mx> in 1 line
                    pos ++;
                    strm >> mnx >> mxx;
                    strm.clear();
                    strm << this -> comlist[pos]; //we hope that file has structure <mn> <mx> in 1 line
                    pos ++;
                    strm >> mny >> mxy;
                    strm.clear();
                    unique_ptr <Control> work (new Control());
                    work -> GenGroup(this -> sz, mnx, mxx, mny, mxy, 0); //base group
                    fstream moved;
                    moved.open("moveX.txt", ios::out | ios::trunc);
                    if(!moved.is_open()){
                        Error(3);
                    }
                    work -> FileGroup();
                    strm << this -> comlist[pos]; //we hope that file has structure <mn> <mx> in 1 line
                    pos ++;
                    strm >> delx;
                    strm.clear();
                    work -> MoveX(delx);
                    gr = work -> RetGroup();
                    for(int i = 0; i < this -> sz; i ++){
                        moved << gr[i].GetX() << " ";
                        moved << gr[i].GetY() << "\n";
                    }
                    moved.close();
                    trg = 0;
                }
                else if(command == "moveY"){//moving Y (DEMO)
                    double dely, mnx, mny, mxx, mxy;
                    vector <Point> gr;
                    strm << this -> comlist[pos]; //we hope that file has structure <mn> <mx> in 1 line
                    pos ++;
                    strm >> mnx >> mxx;
                    strm.clear();
                    strm << this -> comlist[pos]; //we hope that file has structure <mn> <mx> in 1 line
                    pos ++;
                    strm >> mny >> mxy;
                    strm.clear();
                    unique_ptr <Control> work (new Control());
                    work -> GenGroup(this -> sz, mnx, mxx, mny, mxy, 0); //base group
                    fstream moved;
                    moved.open("moveY.txt", ios::out | ios::trunc);
                    if(!moved.is_open()){
                        Error(3);
                    }
                    work -> FileGroup();
                    strm << this -> comlist[pos]; //we hope that file has structure <mn> <mx> in 1 line
                    pos ++;
                    strm >> dely;
                    strm.clear();
                    work -> MoveY(dely);
                    gr = work ->RetGroup();
                    for(int i = 0; i < this -> sz; i ++){
                        moved << gr[i].GetX() << " ";
                        moved << gr[i].GetY() << "\n";
                    }
                    moved.close();
                    trg = 0;
                }
                else if(command == "genfield"){//field creating
                    int numofgroups = 0, tr;
                    double mnx, mny, mxx, mxy;
                    string com = this -> comlist[pos];
                    pos ++;
                    unique_ptr <Field> MyField;
                    unique_ptr <Cluster> MyCluster(new Cluster());
                    unique_ptr <Group> temp;
                    if((com != "get") && (com!="work")){
                        Error(5);//incorrect value
                    } else {
                        if(com == "work"){
                            strm << this -> comlist[pos]; //we hope that file has structure <mn> <mx> in 1 line
                            pos ++;
                            strm >> numofgroups;
                            strm.clear();
                            if(numofgroups < 1){
                                Error(5);
                            }
                            MyField.reset(new Field(numofgroups));
                            for(int i = 0; i < numofgroups; i ++){
                                strm << this -> comlist[pos]; //we hope that file has structure <mn> <mx> in 1 line
                                pos ++;
                                strm >> mnx >> mxx;
                                strm.clear();
                                strm << this -> comlist[pos]; //we hope that file has structure <mn> <mx> in 1 line
                                pos ++;
                                strm >> mny >> mxy;
                                strm.clear();
                                temp.reset(new Group(this -> sz, mnx, mny, mxx, mxy, i));
                                com = "start";
                                while(com != "esc"){
                                    tr = -1;
                                    com = this -> comlist[pos];
                                    pos ++;
                                    if(com == "RN"){//turning of the group (0;0)
                                        double alf;
                                        strm << this -> comlist[pos]; //we hope that file has structure <mn> <mx> in 1 line
                                        pos ++;
                                        strm >> alf;
                                        strm.clear();
                                        alf = GradToRad(alf);
                                        unique_ptr <Control> work = temp -> Ret_Field_Group();
                                        work -> turnNULL(alf);
                                        temp -> Regrupp(work);
                                        tr = 0;
                                    }
                                    else if(com == "RC"){//turning of the group (Center)
                                        double alf;
                                        strm << this -> comlist[pos]; //we hope that file has structure <mn> <mx> in 1 line
                                        pos ++;
                                        strm >> alf;
                                        strm.clear();
                                        alf = GradToRad(alf);
                                        unique_ptr <Control> work = temp -> Ret_Field_Group();
                                        work -> turnCenter(alf);
                                        temp -> Regrupp(work);
                                        tr = 0;
                                    }
                                    else if(com == "MX"){//moving X
                                        double dx;
                                        unique_ptr <Control> work = temp -> Ret_Field_Group();
                                        strm << this -> comlist[pos]; //we hope that file has structure <mn> <mx> in 1 line
                                        pos ++;
                                        strm >> dx;
                                        strm.clear();
                                        work -> MoveX(dx);
                                        temp -> Regrupp(work);
                                        tr = 0;
                                    }
                                    else if(com == "MY"){//moving Y
                                        double dy;
                                        unique_ptr <Control> work = temp -> Ret_Field_Group();
                                        strm << this -> comlist[pos]; //we hope that file has structure <mn> <mx> in 1 line
                                        pos ++;
                                        strm >> dy;
                                        strm.clear();
                                        work -> MoveY(dy);
                                        temp -> Regrupp(work);
                                        tr = 0;
                                    }
                                    if(tr == -1){
                                        if(com != "esc"){
                                            Error(5); //incorrect value
                                        }
                                    }
                                }
                                MyField -> AddGroup(temp);
                            }
                        }
                        else if (com == "get"){
                            MyCluster -> GetFromFile();
                            MyCluster -> ToTxtCluster("Field.txt");
                        }
                        if(com != "get"){
                            MyField -> MakeAllKoord();
                            string ans = this -> comlist[pos];
                            pos ++;
                            MyField -> ToTxt(ans);
                            MyCluster.reset(new Cluster(MyField.get()));
                        }
                        string q;
                        q = this -> comlist[pos];
                        pos ++;
                        if((q != "yes") && (q != "no")){
                            Error(5);
                        }
                        if(q == "no"){
                            q = "esc";
                        }
                        while(q != "esc"){
                            cout << endl << "To Find by Wave enter 'wave'" << endl << "To Find by K-means enter 'km'" << endl
                            << "To Spainning Tree enter 'sptr'" << endl << "To Hierarchy enter 'ie'" << endl
                            << "To Forel alghorithm enter 'fish'" << endl << "To DBSCAN enter 'dbs'" << endl
                            << "To Exit enter 'esc'" << endl;
                            cin >> q;
                            if(q == "wave"){
                                MyCluster -> FindByWaveAlgorithm();
                            } else if (q == "km"){
                                MyCluster -> FindByKMeansAlgorithm();
                            } else if (q == "sptr"){
                                MyCluster -> RunSpainningTreeAlgorithm();
                            } else if (q == "ie"){
                                MyCluster -> FindByHierarchyAlgorithm();
                            } else if (q == "fish"){
                                MyCluster -> FindByForelAlghorithm();
                            } else if (q == "dbs"){
                                MyCluster -> FindByDBSCANAlghorithm();
                            } else if (q != "esc"){
                                Error(5);
                            }
                        }
                        cout <<  endl << this -> comlist[pos];
                        command = this -> comlist[pos];
                        trg = 0;
                    }
                    if(trg == -1){
                        if(command == "exit"){
                            break;
                        } else{
                            Error(5);
                        }
                    }
                    else if(trg == 0){
                        command = this -> comlist[pos];
                        pos ++;
                        trg = -1;
                    }
                }
            }
        }
};

int main(){
    srand(time(NULL));
    setlocale(LC_ALL, "");
    Mylog.open("log.log", ios::out);
    string var = "zzz";
    while((var != "file") && (var != "hand")){
        cout << endl << "Please choose mode ('file'/'hand')" << endl;
        getline(cin, var);
        if((var != "file") && (var != "hand")){
            Error(5);
        }
    }
    if(var == "hand"){
        unique_ptr <Interface> MyProject (new Interface());
        MyProject -> run();
    } else {
        unique_ptr <InterfaceSTR> MyProject (new InterfaceSTR());
        MyProject -> run();
    }
    Mylog.close();
    return 0;
}
