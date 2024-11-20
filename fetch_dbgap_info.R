library(xml2)
library(dplyr)
library(readr)

# if a release is in progress but not completed, this will return the
# in-progress version rather than the released version
fetch_study_xml <- function(this.study) {
    phs.num <- sub("^phs", "", this.study)
    url <- paste0("https://dbgap.ncbi.nlm.nih.gov/ss/dbgapssws.cgi?request=Study&phs=", phs.num)
    read_xml(url)
}

find_study_version <- function(study.xml, version=NULL) {
    if (is.null(version)) {
        return(xml_find_first(study.xml, "Study"))
    } else {
        dat <- study.xml %>%
            xml_find_all("Study")
        for (i in seq_along(dat)) {
            if (xml_attr(dat[i], "v") == version) {
                return(dat[i])
            }
        }
    }
}

xml_version <- function(study.xml, version=NULL) {
    dat <- study.xml %>%
        find_study_version(version)
    tibble(
        version=xml_attr(dat, "v"),
        participant_set=xml_attr(dat, "p")
    )
}

xml_release_date <- function(study.xml, version=NULL) {
    stat <- study.xml %>%
        find_study_version(version) %>%
        xml_find_first("Status")
    substr(xml_attr(stat, "release_date"), 1, 10)
}

xml_consent <- function(study.xml, version=NULL) {
    cons <- study.xml %>%
        find_study_version(version) %>%
        xml_find_first("Policy") %>%
        xml_find_all("ConsentGroup")
    tibble(consent_code=xml_attr(cons, "name"),
           consent_group=xml_attr(cons, "title"),
           DUL=xml_text(xml_child(cons, "Use-Restriction")),
           DAC=xml_attr(cons, "dac_name"))
}

xml_acknowledgements <- function(study.xml, version=NULL) {
    phs <- study.xml %>%
        xml_find_first("Study") %>%
        xml_attr("phs")
    pol <- study.xml %>%
        find_study_version(version) %>%
        xml_find_first("Policy")
    ack <- xml_text(xml_child(pol, "acknowledgement_statement"))
    if (ack == "") {
        ack <- xml_child(pol, "DUC_AcknowledgementStatement")
        if (length(xml_child(ack, "File")) > 0) {
            vers <- xml_version(study.xml, version)
            ack <- sprintf("https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/GetAcknowledgementStatement.cgi?study_id=phs%s.v%s.p%s", phs, vers$version, vers$participant_set)
        }
    }
    return(ack)
}

# Determine if a study version is marked as "sensitive" for GSR posting.
xml_gsr <- function(study.xml, version=NULL) {
    gsr <- study.xml %>%
        find_study_version(version) %>%
        xml_find_first("Policy") %>%
        xml_find_all("GSR_Access")
    xml_attr(gsr, "gsr_mode_label")
}

xml_partners <- function(study.xml, version=NULL) {
    part <- study.xml %>%
        find_study_version(version) %>%
        xml_find_first("Policy") %>%
        xml_find_first("TrustedPartners") %>%
        xml_find_all("TrustedPartner")
    tibble(trusted_partner=xml_attr(part, "trp_db_name"))
}

xml_info <- function(study.xml, version=NULL) {
    dat <- study.xml %>%
        find_study_version(version)
    info <- dat %>%
        xml_find_first("StudyInfo")
    tibble(
        dbgap_abbrev=xml_attr(dat, "handle"),
        dbgap_name=xml_text(xml_find_all(info, "StudyNameEntrez")),
        accession=xml_attr(info, "accession"),
        parent=xml_attr(info, "parentAccession"),
        n_participants=xml_attr(dat, "num_participants")
    ) %>%
        mutate(parent=ifelse(accession == parent, NA, parent))
}

xml_children <- function(study.xml, version=NULL) {
    dat <- study.xml %>%
        find_study_version(version)
    info <- dat %>%
        xml_find_first("StudyInfo")
    sapply(xml_find_all(info, "childAccession"), xml_text)
}

xml_pi <- function(study.xml, version=NULL) {
    people <- study.xml %>%
        find_study_version(version) %>%
        xml_find_first("Authority") %>%
        xml_find_first("Persons") %>%
        xml_find_all("Person")
    roles <- xml_find_all(people, "Role") %>%
        xml_text()
    pi <- people[which(roles == "PI")]
    tibble(PI_name=paste(xml_attr(pi, "fname"), xml_attr(pi, "lname")),
           PI_email=xml_attr(pi, "email"))
}

xml_data <- function(study.xml, version=NULL) {
    dat <- study.xml %>%
        find_study_version(version) %>%
        xml_find_first("StudyInfo") %>%
        xml_find_first("StudyTypes")
    tibble(data_types=xml_text(xml_children(dat)))
}


parse_dbgap_consent <- function(study.xml) {
    lapply(names(study.xml), function(x) {
        tryCatch({
            xml_consent(study.xml[[x]]) %>%
                mutate(phs=x) %>%
                filter(!grepl("^Exchange Area", consent_group))
        }, error=function(e) tibble(phs=x)
        )
    }) %>%
    bind_rows()
}


parse_dbgap_info <- function(study.xml) {
    lapply(names(study.xml), function(x) {
        tryCatch({
            xml_info(study.xml[[x]]) %>%
                mutate(phs=x)
        }, error=function(e) tibble(phs=x)
        )
    }) %>%
        bind_rows()
}


parse_dbgap_children <- function(study.xml) {
    lapply(names(study.xml), function(x) {
        tryCatch({
            xml_children(study.xml[[x]]) %>%
                substr(1,9)
        }, error=function(e) NULL
        )
    }) %>%
        setNames(names(study.xml))
}


parse_dbgap_pi <- function(study.xml) {
    lapply(names(study.xml), function(x) {
        tryCatch({
            xml_pi(study.xml[[x]]) %>%
                mutate(phs=x)
        }, error=function(e) tibble(phs=x)
        )
    }) %>%
        bind_rows()
}


parse_dbgap_partners <- function(study.xml) {
    lapply(names(study.xml), function(x) {
        tryCatch({
            xml_partners(study.xml[[x]]) %>%
                mutate(phs=x)
        }, error=function(e) tibble(phs=x)
        )
    }) %>%
        bind_rows()
}


fetch_telemetry <- function(this.study) {
    url <- paste0("https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/GetSampleStatus.cgi?study_id=", this.study, "&rettype=xml")
    read_xml(url)
}

count_samples <- function(study.xml) {
    study.xml %>%
        xml_find_first("Study") %>%
        xml_find_first("SampleList") %>%
        xml_length()
}

count_samples_html <- function(this.study) {
    url <- paste0("https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/GetSampleStatus.cgi?study_id=", this.study, "&rettype=html")
    html <- scan(url, what="", sep="\n", nmax=50, quiet=TRUE)
    html[which(grepl("sample ids are Loaded", html))] %>%
        stringr::str_extract("[:digit:]+") %>%
        as.integer()
}


fetch_dars <- function(this.study, v, p) {
    url <- paste0("https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/GetAuthorizedRequestDownload.cgi?study_id=", this.study, ".v", v, ".p", p)
    dars <- suppressWarnings(read_tsv(url, show_col_types = FALSE))
    if (nrow(dars) == 0) dars <- tibble()
    dars
}

count_dars <- function(this.study, v, p) {
    url <- paste0("https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/GetAuthorizedRequestDownload.cgi?study_id=", this.study, ".v", v, ".p", p)
    dars <- suppressWarnings(readLines(url))
    num <- if (length(dars) > 1) (length(dars) - 1) else 0
    return(num)
}
