import argparse

def inputHandling():
    parser = argparse.ArgumentParser(description='Gets the SV consensus.')

    parser.add_argument('-f', '--sv_folder', help='Folder containing folders of samples with raw outputs from SV callers (comma-separated). More information on the structure of the samples folder in readme.', required=True)

    parser.add_argument('-mod', '--model', help='Model used for SV discovery (default default.model).', required=False, default="default.model")

    parser.add_argument('-o', '--output', help='Output file prefix.', default="consensuSV_")

    parser.add_argument('-m', '--min_overlap', help='Minimum number of SVs in the neighbourhood for the SV to be reported (default 3).', type=int, default=3)

    parser.add_argument('-t', '--train', help='Creates new model. Requires truth.vcf to be present in all the sv folders. VCF file truth.vcf is preprocessed even if flag --no_preprocess is set.', action="store_true", required=False)

    parser.add_argument('-np', '--no_preprocess', help='Flag used for skipping the preprocessing process - all the preprocessed files should be in temp/ folder.', action="store_true", required=False)

    args = parser.parse_args()

    return args