# ConsensuSV
A tool for getting consensus of SVs.

Options:

Short option | Long option | Description
-------------- | --------------- | ---------------
-f | --sv_folder | Folder containing folders of samples with raw outputs from SV callers (comma-separated). More information on the structure of the samples folder in next paragraph.
-mod | --model | Model used for SV discovery (default: default.model).
-o | --output | Output file prefix.
-m | --min_overlap | Minimum number of SVs in the neighbourhood for the SV to be reported (default: 3).
-t | --train | Creates new model. Requires truth.vcf to be present in all the sv folders. VCF file truth.vcf is preprocessed even if flag --no_preprocess is set. If the model is trained, it is required to rerun the program to get the consensus.
-np | --no_preprocess | Flag used for skipping the preprocessing process - all the preprocessed   files should be in temp/ folder.
