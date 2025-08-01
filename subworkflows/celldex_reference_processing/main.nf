include { CELLTYPES_CELLDEXDOWNLOAD } from '../../modules/celltypes/celldexdownload'

workflow CELLDEX_REFERENCE_PROCESSING {
    take:
    reference_string

    main:
    def reftars = Channel.empty()
    def ref_list = reference_string.map{it -> it.toString().split("--")}.flatten()
    log.info "Processing references: ${ref_list}"
    def to_download = []
    // For each reference in the ref_list channel, check if it ends in .tar.gz, if so, add it to the reftars channel.
    // If it does not end in .tar.gz, add it to the to_download list
    ref_list.each { ref ->
        if (ref.endsWith(".tar.gz")) {
            reftars = reftars.mix(file(ref))
        } else if (ref.endsWith(".tar")) {
            reftars = reftars.mix(file(ref + ".tar.gz"))
        } else {
            to_download.add(ref)
        }
    }

    if (to_download.size() > 0) {
        Channel.fromList(to_download) | CELLTYPES_CELLDEXDOWNLOAD
        reftars = reftars.mix(CELLTYPES_CELLDEXDOWNLOAD.out.tar)
    }
    emit: referenceTars = reftars.collect()
}
