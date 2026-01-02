#pragma once

#include "streamlineAutoencoder.h"

namespace NIBR
{

    // inp:   path to input file containing latent representations of streamlines
    // out:   path to output tractogram file to save
    // force: if true, overwrites out if it already exists
    bool decodeAndSave(std::string inp, std::string out, bool force, StreamlineAutoencoder& model, int batchSize, float resampleStepSize = 0);

    // inp:   path to input tractogram file
    // out:   path to output file containing latent representations of streamlines
    // force: if true, overwrites out if it already exists
    bool encodeAndSave(std::string inp, std::string out, bool force, StreamlineAutoencoder& model, int batchSize, bool performResampling = false);

}