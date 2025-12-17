#include "dMRI/tractography/tractogram.h"
#include "StreamlineAutoencoder.h"
using namespace NIBR;

bool encodeAndSave(std::string inp, std::string out, bool force, StreamlineAutoencoder& model, int batchSize);