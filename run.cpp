#include "Count.h"

using namespace std;

int main(int argc, char* argv[]){
	init_workers();
	SetLCR_Construction(  "/blogel/test.graph",
				          "/blogel/p2",
			              "/blogel/index_",
			        1049, 123980, 32);

	worker_finalize();
	

	return 0;
}
