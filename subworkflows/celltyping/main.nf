include { CELLTYPES_SINGLER } from '../../modules/celltypes/singler/'
workflow CELLTYPE_ASSIGNMENT {
    take:
    sampleName
    h5ad
    celldex_ref_dirs
    celldex_reference_label
    
    main:

    CELLTYPES_SINGLER(sampleName, h5ad, celldex_ref_dirs, celldex_reference_label)

}
