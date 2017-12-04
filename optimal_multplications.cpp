#include <iostream>
#include <vector>
#include <limits>
#include <fstream>

using namespace std;

/**
* Compute optimal ordering of matrix multiplication.
* c contains the number of columns for each of the n matrices.
* c[ 0 ] is the number of rows in matrix 1.
* The minimum number of multiplications is left in m[ 1 ][ n ].
* Actual ordering is computed via another procedure using lastChange.
* m and lastChange are indexed starting at 1, instead of 0.
* Note: Entries below main diagonals of m and lastChange
* are meaningless and uninitialized.
*/
//void optMatrix( const vector<int> & c, matrix<int> & m, matrix<int> & lastChange ) {
void optMatrix( const vector<int> & c, vector< vector<long> > & m, vector< vector<int> > & lastChange ) {
  int n = c.size() - 1;

  for( int left = 1; left <= n; ++left )
    m[ left ][ left ] = 0;

  for( int k = 1; k < n; ++k ) {  // k is right - left
    for( int left = 1; left <= n - k; ++left ) {
      // For each position
      int right = left + k;
      m[ left ][ right ] = std::numeric_limits<long>::max();  // set to "infinity"
      for( int i = left; i < right; ++i ) {
        long thisCost = m[ left ][ i ] + m[ i + 1 ][ right ] + c[ left - 1 ] * c[ i ] * c[ right ];
//cout<< /*left << " x " << i+1*/ "m"<<left<<"," <<right << " = " <<thisCost << "   From: m[" << left << "][" << i << "] + m["<<i+1<<"]["<<right<<"] + c["<<left-1<<"] * c["<<i<<"] * c["<<right<<"]   left: "<<left<<" right: "<<right <<endl;
        if( thisCost < m[ left ][ right ] ) {  // Update min
          m[ left ][ right ] = thisCost;
          lastChange[ left ][ right ] = i;
        }
      }  
//cout<<"==============="<<endl;
    }
  }


  return;
}

bool extractContentFromLine( const string & line, int & value, const int & line_num ) {
  if( line == "" ) {
    cout << "Format Error at line " << line_num << ": Empty line" << endl;
    return false;
  }
   // Checks line if any non number character, if so, terminate input process.
  for( std::string::const_iterator itr = line.begin(); itr != line.end(); ++itr ) {
    if( *itr < '0' or *itr > '9' ) {
      cout << "Format Error at line " << line_num << ": Non-number character" << endl;
      return false;
    }
  }
  
  value = std::stoi( line );
  return true;
}

 // Takes in string file name and a vector, will read matrix formated text file.
bool readFile( const string & filename, vector<int> & dimensions ) {
  ifstream file( filename.c_str() );
  if( file.fail() ) {
    cout << "Unable to open file..." << endl;
    return false;
  } 
  
  string line = "";

  int counter = 1;
  int value = 0;
  while( std::getline( file, line ) ) {
    //cout << line << endl;
    if( !extractContentFromLine( line, value, counter )) {
       return false;
    }
    ++counter;
//cout<<"val:"<<value<<endl;
    dimensions.push_back(value);
  }
  
  return true;
}

void printOptMatrix( const vector< vector<long> > & m ) {
  cout << "  ";
  for( unsigned int i = 1; i < m.size(); ++i )
    cout << i << " ";
  cout << endl;
  cout << " -----------------------------------------" << endl;  

  for( unsigned int i = 1; i < m.size(); ++i ) {
    cout << i << "|";
    for( unsigned int j = 1; j < m[i].size(); ++j )
      cout << m[i][j] << " ";
    cout << endl;
  }
  cout << " -----------------------------------------" << endl;  
}

void printOptMatrixLastChange( const vector< vector<int> > & lastChange ) {
  cout << "  ";
  for( unsigned int i = 1; i < lastChange.size(); ++i )
    cout << i << " ";
  cout << endl;
  cout << " -----------------------------------------" << endl;  

  for( unsigned int i = 1; i < lastChange.size(); ++i ) {
    cout << i << "|";
    for( unsigned int j = 1; j < lastChange[i].size(); ++j )
      cout << lastChange[i][j] << " ";
    cout << endl;
  }
  cout << " -----------------------------------------" << endl;  
}

 // helper to below fct
void traverseWithinMatrixMult( const vector<int> & dimensions, const vector< vector<int> > & lastChange, const int & left, const int & right, string & equation ) {
  if( left == right ) {
    //equation += "(";
    //equation += to_string( left );
    equation += to_string( dimensions[left] );
    //equation += ")";
    return;
  }
  
  int cutoff = lastChange[left][right];  //( left...cutoff) ( cuttoff+1 ....right)

  equation += "(";
  traverseWithinMatrixMult( dimensions, lastChange, left, cutoff, equation );
  if( right - left == 1 )
    equation += "*";
  traverseWithinMatrixMult( dimensions, lastChange, cutoff+1, right, equation );
  equation += ")";
  
  return;
}

 // recursive solution to print out best order of multiplication. M(left,right)
void printBestOrder( const vector<int> & dimensions, const vector< vector<int> > & lastChange, const int & left, const int & right ) {
  string equation = "";
  traverseWithinMatrixMult( dimensions, lastChange, left, right, equation );
  cout << equation << endl;

  return; 
}

int main( int argc, char **argv ) {
  if( argc != 2 ) {
    cout << "Usage: ./optimal_multplications [Matrix Text File, typically dimensions_file,txt]" << endl;
    return 0;
  }

  const string dimension_file( argv[1] );
  vector<int> dimensions;

  if( !readFile( dimension_file, dimensions ) ) {
    cout << "Failed to input file..." << endl;
    return 0;
  }

  if( dimensions.size() == 0 ) {
    cout << "Nothing was in the file for the program to populate..." << endl;
    return 0;
  }

  vector< vector<long> > m( dimensions.size(), vector<long>( dimensions.size(), 0 ) );
  vector< vector<int> > last_change( dimensions.size(), vector<int>( dimensions.size(), 0 ) );

  optMatrix( dimensions, m, last_change );
  //printOptMatrix( m );
  //printOptMatrixLastChange( last_change );

  cout << "Optimized number of multiplications: " << m[ 1 ][ m.size() - 1 ] << endl;
  cout << "Correct order of multiplication to achieve the optimized number: ";
  printBestOrder( dimensions ,last_change, 1, dimensions.size()-1 );
  cout << endl;
  cout << "NOTE: The final multiplied result from the correct order does neccessarily equal the optimized number of multiplications. The optimized number of multiplications simply means that it takes this many steps to get the result of multiplying N number of Matrices." << endl;

  return 0;
}
