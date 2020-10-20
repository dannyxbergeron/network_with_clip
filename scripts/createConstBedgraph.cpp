#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

using namespace std;


void getCoord(const string& filename, vector<string>& coord)
{
    ifstream infile(filename);
    string line;
    int start(0), end, tmpStart, tmpEnd;

    while(getline(infile, line))
    {
        for(int i = 0; i < 3; i++)
        {
            end = line.find('\t', start);
            start = end + 1;
        }
        tmpEnd = end;

        for(int i = 0; i < 2; i++)
        {
            tmpStart = tmpEnd + 1;
            tmpEnd = line.find('\t', tmpStart);
        }

        if(line.substr(tmpStart,tmpEnd - tmpStart) == "transcript")
        {
            coord.push_back(line.substr(0, end));
        }
        start = 0;
    }
}

bool splitString(const string& s, string& chr, int& start, int& end)
{
    string col;
    int beg(0), find;
    find = s.find('\t');
    chr = s.substr(0, find);

    beg = find + 1;
    find = s.find('\t', beg);
    start = stoi(s.substr(beg, find - beg));

    beg = find + 1;
    find = s.find('\t', beg);
    end = stoi(s.substr(beg, find - beg));
}


void keepTranscriptInfo(vector<string>& coord)
{
    string chr, tChr("chr0");
    int start, end, tStart(-1), tEnd(-1);
    float val;
    int pos = -1;

    vector<string> chrNames = {
        "chr0", "chr1", "chr10", "chr11", "chr12", "chr13","chr14", "chr15",
        "chr16", "chr17", "chr18", "chr19", "chr2", "chr20", "chr21", "chr22",
        "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chrX", "chrY"
    };
    unordered_map<string, int> chrs;
    for(int i = 0; i < chrNames.size(); i++)
    {
        chrs[chrNames[i]] = i;
    }


    ios_base::sync_with_stdio(false);

    while(cin >> chr >> start >> end >> val)
    {
        if(chrs[tChr] > chrs[chr])
            continue;

        if(start > tEnd || chrs[tChr] < chrs[chr])
        {
            pos++;
            if(pos == coord.size())
                return;
            splitString(coord[pos], tChr, tStart, tEnd);
        }

        if(chr == tChr && start > tStart && end < tEnd)
        {
            cout << chr << '\t' << start << '\t' << end << '\t' << val << endl;
        }
    }
}


int main(int argc, char **argv)
{
    vector<string> coord;
    getCoord(argv[1], coord);


    keepTranscriptInfo(coord);


    return 0;
}
