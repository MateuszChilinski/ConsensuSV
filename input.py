import argparse

def inputHandling():
    parser = argparse.ArgumentParser(description='Gets the SV consensus.')
    parser.add_argument('-f', '--sv_folder', help='Folder containing raw outputs from SV callers.', required=True)

    parser.add_argument('-s', '--sample', help='Name of the sample.', required=True)

    parser.add_argument('-o', '--output', help='Output file.', default="output.vcf")

    parser.add_argument('-m', '--min_overlap', help='Minimum number of SVs in the neighbourhood for the SV to be reported (default 3).', type=int, default=3)

    parser.add_argument('-t', '--truth', help='File used for training new model.', required=False)

    parser.add_argument('-np', '--no_preprocess', help='Flag used for skipping the preprocessing process - all the preprocessed files should be in temp/ folder.', action="store_true", required=False)

    args = parser.parse_args()

    return args