import os

import sjcab_peak2anno_db as db
import sjcab_peak2anno_db._chromhmm as chromhmm
import sjcab_peak2anno_db._cli as cli
import sjcab_peak2anno_db._dedup as dedup
import sjcab_peak2anno_db._download as download
import sjcab_peak2anno_db._gencode as gencode
import sjcab_peak2anno_db._install as install
import sjcab_peak2anno_db._regions as regions
from sjcab_peak2anno_db._gencode import ensembl_gtf_url
import sjcab_peak2anno_db._segway as segway


def _write_mini_gencode_region_gtf(tmp_path):
    gtf = tmp_path / "mini.gtf"
    gtf.write_text(
        "\n".join(
            [
                "##sequence-region chr1 1 600",
                'chr1\tGENCODE\tgene\t101\t220\t.\t+\t.\tgene_id "ENSG1.2"; gene_name "GENE1"; gene_type "protein_coding";',
                'chr1\tGENCODE\ttranscript\t101\t220\t.\t+\t.\tgene_id "ENSG1.2"; transcript_id "ENST1.3"; gene_name "GENE1"; gene_type "protein_coding";',
                'chr1\tGENCODE\texon\t101\t120\t.\t+\t.\tgene_id "ENSG1.2"; transcript_id "ENST1.3"; gene_name "GENE1"; gene_type "protein_coding";',
                'chr1\tGENCODE\texon\t201\t220\t.\t+\t.\tgene_id "ENSG1.2"; transcript_id "ENST1.3"; gene_name "GENE1"; gene_type "protein_coding";',
                'chr1\tGENCODE\ttranscript\t101\t220\t.\t+\t.\tgene_id "ENSG1.2"; transcript_id "ENST1.4"; gene_name "GENE1"; gene_type "protein_coding";',
                'chr1\tGENCODE\texon\t101\t150\t.\t+\t.\tgene_id "ENSG1.2"; transcript_id "ENST1.4"; gene_name "GENE1"; gene_type "protein_coding";',
                'chr1\tGENCODE\texon\t181\t220\t.\t+\t.\tgene_id "ENSG1.2"; transcript_id "ENST1.4"; gene_name "GENE1"; gene_type "protein_coding";',
                'chr1\tGENCODE\tgene\t401\t500\t.\t-\t.\tgene_id "ENSG2.1"; gene_name "GENE2"; gene_type "protein_coding";',
                'chr1\tGENCODE\ttranscript\t401\t500\t.\t-\t.\tgene_id "ENSG2.1"; transcript_id "ENST2.1"; gene_name "GENE2"; gene_type "protein_coding";',
                'chr1\tGENCODE\texon\t401\t450\t.\t-\t.\tgene_id "ENSG2.1"; transcript_id "ENST2.1"; gene_name "GENE2"; gene_type "protein_coding";',
                'chr1\tGENCODE\texon\t481\t500\t.\t-\t.\tgene_id "ENSG2.1"; transcript_id "ENST2.1"; gene_name "GENE2"; gene_type "protein_coding";',
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return gtf


def _chromhmm_data_js():
    return """
var data_epg=
[
  {
    "eid":"E003",
    "group":"Blood & T-cell",
    "mnemonic":"BLD.CD19",
    "name":"GM12878 Lymphoblastoid Cells",
    "name_EDACC_R9":"GM12878/Lymphoblastoid Cell-Line",
    "class":"Class 1",
    "info":""
  },
  {
    "eid":"E063",
    "group":"Brain / Fetal",
    "mnemonic":"BRN.FETAL",
    "name":"Fetal Brain Male",
    "name_EDACC_R9":"Fetal Brain Male",
    "class":"Class 2",
    "info":""
  }
];
"""


def _segway_index_html():
    return """
<a href="interpreted">Directory</a>
<a href=segway_encyclopedia.bed.gz>BED format</a>
<a href=label_info.tab>tab-delimited format</a>
<a href="caas.bed.gz">BED format</a>
"""


def _segway_interpreted_html():
    return """
<a href="GM12878.bed.gz">GM12878.bed.gz</a>
<a href="H1-HESC.bed.gz">H1-HESC.bed.gz</a>
<a href="BRAIN_ANGULAR_GYRUS.bed.gz">BRAIN_ANGULAR_GYRUS.bed.gz</a>
"""


def _segway_encode_file(accession, href, submitted_file_name, **kwargs):
    payload = {
        "@type": ["File", "Item"],
        "accession": accession,
        "href": href,
        "file_format": "bed",
        "file_format_type": "bed9",
        "output_type": "semi-automated genome annotation",
        "assembly": "GRCh37",
        "status": "released",
        "submitted_file_name": submitted_file_name,
        "aliases": ["encode:segway-recolored"],
    }
    payload.update(kwargs)
    return payload


def test_download_file_reports_every_five_percent(monkeypatch, tmp_path):
    class FakeResponse:
        headers = {"Content-Length": "100"}

        def __init__(self):
            self.remaining = b"x" * 100

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def read(self, size):
            chunk = self.remaining[:size]
            self.remaining = self.remaining[size:]
            return chunk

    def fake_urlopen(url, timeout):
        assert url.full_url == "https://example.org/file.bed"
        assert url.headers["User-agent"] == "sjcab-peak2anno-db"
        assert timeout == 120
        return FakeResponse()

    monkeypatch.setattr(download.urllib.request, "urlopen", fake_urlopen)
    messages = []

    download.download_file(
        "https://example.org/file.bed",
        tmp_path / "file.bed",
        progress=messages.append,
    )

    assert messages[0] == "download file.bed:"
    assert messages[1:-1] == [
        "{}..".format(percent) for percent in range(5, 105, 5)
    ]
    assert messages[-1] == "done"
    assert (tmp_path / "file.bed").read_text(encoding="utf-8") == "x" * 100


def test_cli_progress_appends_download_percent_fragments(capsys):
    cli._stderr_progress("download E063_18_core_K27ac_dense.bed.gz:")
    cli._stderr_progress("5..")
    cli._stderr_progress("10..")
    cli._stderr_progress("done")

    captured = capsys.readouterr()
    assert captured.err == "download E063_18_core_K27ac_dense.bed.gz:5..10..done\n"


def test_default_version_uses_registry_defined_version_not_largest():
    assert db.versions("hg38", "gene") == ("v31",)
    assert db.default_version("hg38", "gene") == "v31"
    assert db.resource("hg38", "gene").name == "gencode.v31.hg38.gene.bed.withtype"

    assert db.versions("hg19", "gene") == ("v31lift37",)
    assert db.default_version("hg19", "gene") == "v31lift37"
    assert (
        db.resource("hg19", "gene").name
        == "gencode.v31lift37.hg19.gene.bed.withtype"
    )


def test_legacy_deduplong_annotation_resolves_to_isoform_set(tmp_path):
    target = db.installed_path("hg38", "deduplong", data_dir=tmp_path)

    assert target == tmp_path / "bed" / "hg38" / "def" / "deduplong.gene.bed"


def test_refresh_default_links_to_registry_defined_version(tmp_path):
    entry = next(
        resource
        for resource in db.iter_resources()
        if resource.species == "hg38"
        and resource.isoform_set == "all"
        and resource.annotation == "gene"
        and resource.version == db.default_version("hg38")
    )
    version_path = db.installed_path(
        entry.species,
        entry.annotation,
        entry.version,
        data_dir=tmp_path,
        isoform_set=entry.isoform_set,
    )
    version_path.parent.mkdir(parents=True, exist_ok=True)
    version_path.write_text("defined default\n", encoding="utf-8")

    install._refresh_default(entry, tmp_path)

    default_path = db.installed_path(
        "hg38", "gene", "default", data_dir=tmp_path, isoform_set="all"
    )
    assert default_path.read_text(encoding="utf-8") == "defined default\n"
    default_dir = tmp_path / "bed" / "hg38" / "def"
    if default_dir.is_symlink():
        assert os.readlink(str(default_dir)) == "v31"


def test_install_gencode_bed_uses_bed_component_layout(tmp_path):
    entry = next(
        resource
        for resource in db.iter_resources()
        if resource.species == "sacCer3"
        and resource.version == db.default_version("sacCer3")
        and resource.isoform_set == "all"
        and resource.annotation == "gene"
    )

    install._install_gene(entry, tmp_path, overwrite=True)
    install._refresh_default(entry, tmp_path)

    version_path = tmp_path / "bed" / "sacCer3" / entry.version / "all.gene.bed"
    default_dir = tmp_path / "bed" / "sacCer3" / "def"
    assert version_path.exists()
    if default_dir.is_symlink():
        assert os.readlink(str(default_dir)) == entry.version
    assert (default_dir / "all.gene.bed").exists()


def test_install_gencode_bed_cli_reuses_existing_files_by_default(
    monkeypatch, tmp_path
):
    called = {}

    def fake_install(data_dir=None, overwrite=True):
        called["data_dir"] = data_dir
        called["overwrite"] = overwrite
        return tmp_path

    monkeypatch.setattr(cli, "install_gencode_beds", fake_install)

    assert cli.main(["install-bed", "-d", str(tmp_path)]) == 0
    assert called == {"data_dir": str(tmp_path), "overwrite": False}


def test_install_gencode_bed_cli_accepts_explicit_overwrite(
    monkeypatch, tmp_path
):
    called = {}

    def fake_install(data_dir=None, overwrite=True):
        called["data_dir"] = data_dir
        called["overwrite"] = overwrite
        return tmp_path

    monkeypatch.setattr(cli, "install_gencode_beds", fake_install)

    assert cli.main(["install-bed", "-d", str(tmp_path), "--overwrite"]) == 0
    assert called == {"data_dir": str(tmp_path), "overwrite": True}


def test_download_gencode_bed_accepts_short_output_dir(monkeypatch, tmp_path):
    called = {}

    def fake_download(species, version, output_dir, **kwargs):
        called["args"] = (species, version, output_dir)
        called["kwargs"] = kwargs
        return tmp_path / "genes.bed"

    monkeypatch.setattr(cli, "download_and_convert_gencode_gtf", fake_download)

    assert (
        cli.main(["download-bed", "hg38", "v31", "-o", str(tmp_path)])
        == 0
    )
    assert called["args"] == ("hg38", "v31", str(tmp_path))


def test_path_install_reuses_existing_files_by_default(
    monkeypatch, tmp_path, capsys
):
    called = {}

    def fake_install(data_dir=None, overwrite=True):
        called["data_dir"] = data_dir
        called["overwrite"] = overwrite
        return tmp_path

    def fake_path(*args, **kwargs):
        return tmp_path / "resolved.bed"

    monkeypatch.setattr(cli, "install_data", fake_install)
    monkeypatch.setattr(cli, "path", fake_path)

    assert cli.main(["path", "hg38", "gene", "-d", str(tmp_path), "-I"]) == 0
    assert called == {"data_dir": str(tmp_path), "overwrite": False}
    captured = capsys.readouterr()
    assert captured.out == "{}\n".format(tmp_path / "resolved.bed")


def test_download_gencode_bed_defaults_to_current_output_dir(monkeypatch):
    called = {}

    def fake_download(species, version, output_dir, **kwargs):
        called["args"] = (species, version, output_dir)
        called["kwargs"] = kwargs
        return "genes.bed"

    monkeypatch.setattr(cli, "download_and_convert_gencode_gtf", fake_download)

    assert cli.main(["download-bed", "hg38", "v31"]) == 0
    assert called["args"] == ("hg38", "v31", ".")


def test_download_gencode_gtf_reuses_existing_output_dir_file(monkeypatch, tmp_path):
    existing = tmp_path / "gencode.v31.annotation.gtf.gz"
    existing.write_text("already downloaded\n", encoding="utf-8")

    def fail_download(url, destination, **kwargs):
        raise AssertionError("download should not be called")

    monkeypatch.setattr(gencode, "_download_file", fail_download)

    output = db.download_gencode_gtf("hg38", "v31", tmp_path, overwrite=True)

    assert output == existing
    assert existing.read_text(encoding="utf-8") == "already downloaded\n"


def test_download_gencode_gtf_reuses_existing_working_dir_file(monkeypatch, tmp_path):
    working_dir = tmp_path / "work"
    output_dir = tmp_path / "out"
    working_dir.mkdir()
    output_dir.mkdir()
    existing = working_dir / "gencode.v31.annotation.gtf.gz"
    existing.write_text("already downloaded\n", encoding="utf-8")

    def fail_download(url, destination, **kwargs):
        raise AssertionError("download should not be called")

    monkeypatch.setattr(gencode, "_download_file", fail_download)
    monkeypatch.chdir(working_dir)

    output = db.download_gencode_gtf("hg38", "v31", output_dir, overwrite=True)

    assert output == existing
    assert existing.read_text(encoding="utf-8") == "already downloaded\n"


def test_download_gencode_bed_writes_version_layout(tmp_path):
    gtf = _write_mini_gencode_region_gtf(tmp_path)

    output = db.download_and_convert_gencode_gtf(
        "hg38",
        "v31",
        tmp_path,
        gtf_path=gtf,
    )

    version_dir = tmp_path / "bed" / "hg38" / "v31"
    assert output == version_dir / "all.gene.bed"
    assert (version_dir / "deduplong.gene.bed").exists()
    assert not (version_dir / "all.tss.bed").exists()
    assert not (version_dir / "all.tes.bed").exists()
    assert not (version_dir / "deduplong.tss.bed").exists()
    assert not (version_dir / "deduplong.tes.bed").exists()
    default_dir = tmp_path / "bed" / "hg38" / "def"
    if default_dir.is_symlink():
        assert os.readlink(str(default_dir)) == "v31"
    assert (default_dir / "all.gene.bed").exists()


def test_install_accepts_component_names(monkeypatch, tmp_path):
    called = {}

    def fake_install(data_dir=None, overwrite=True, components=None, progress=None):
        called["data_dir"] = data_dir
        called["overwrite"] = overwrite
        called["components"] = components
        called["progress"] = progress
        return tmp_path

    monkeypatch.setattr(cli, "install_data", fake_install)

    assert (
        cli.main(
            [
                "install",
                "gencode-bed",
                "cgi",
                "--data-dir",
                str(tmp_path),
                "--no-overwrite",
            ]
        )
        == 0
    )
    assert called == {
        "data_dir": str(tmp_path),
        "overwrite": False,
        "components": ["gencode-bed", "cgi"],
        "progress": cli._stderr_progress,
    }


def test_cli_install_without_components_uses_default_selection(monkeypatch, tmp_path):
    called = {}

    def fake_install(data_dir=None, overwrite=True, components=None, progress=None):
        called["data_dir"] = data_dir
        called["overwrite"] = overwrite
        called["components"] = components
        called["progress"] = progress
        return tmp_path

    monkeypatch.setattr(cli, "install_data", fake_install)

    assert cli.main(["install", "-d", str(tmp_path)]) == 0
    assert called == {
        "data_dir": str(tmp_path),
        "overwrite": False,
        "components": None,
        "progress": cli._stderr_progress,
    }


def test_cli_install_accepts_explicit_overwrite(monkeypatch, tmp_path):
    called = {}

    def fake_install(data_dir=None, overwrite=True, components=None, progress=None):
        called["data_dir"] = data_dir
        called["overwrite"] = overwrite
        called["components"] = components
        called["progress"] = progress
        return tmp_path

    monkeypatch.setattr(cli, "install_data", fake_install)

    assert cli.main(["install", "-d", str(tmp_path), "--overwrite"]) == 0
    assert called == {
        "data_dir": str(tmp_path),
        "overwrite": True,
        "components": None,
        "progress": cli._stderr_progress,
    }


def test_install_accepts_short_component_options(monkeypatch, tmp_path):
    called = {}

    def fake_install(data_dir=None, overwrite=True, components=None, progress=None):
        called["data_dir"] = data_dir
        called["overwrite"] = overwrite
        called["components"] = components
        called["progress"] = progress
        return tmp_path

    monkeypatch.setattr(cli, "install_data", fake_install)

    assert (
        cli.main(
            [
                "install",
                "-c",
                "gencode-bed",
                "-c",
                "cgi",
                "-d",
                str(tmp_path),
                "-n",
            ]
        )
        == 0
    )
    assert called == {
        "data_dir": str(tmp_path),
        "overwrite": False,
        "components": ["gencode-bed", "cgi"],
        "progress": cli._stderr_progress,
    }


def test_install_data_default_installs_expected_components(monkeypatch, tmp_path):
    calls = []

    def fake_gencode_beds(data_dir=None, overwrite=True):
        calls.append(("gencode-bed", data_dir, overwrite))
        return tmp_path

    def fake_gencode_features(data_dir=None, overwrite=True, progress=None):
        calls.append(("gencode-feature", data_dir, overwrite))
        return tmp_path

    def fake_blacklists(data_dir=None, overwrite=True):
        calls.append(("blacklists", data_dir, overwrite))
        return tmp_path

    def fake_cgi(data_dir=None, overwrite=True):
        calls.append(("cgi", data_dir, overwrite))
        return tmp_path

    def fake_manifest(target_root):
        calls.append(("manifest", target_root, None))

    monkeypatch.setattr(install, "install_gencode_beds", fake_gencode_beds)
    monkeypatch.setattr(install, "install_gencode_features", fake_gencode_features)
    monkeypatch.setattr(install, "install_blacklists", fake_blacklists)
    monkeypatch.setattr(install, "install_cgi", fake_cgi)
    monkeypatch.setattr(install, "_write_manifest", fake_manifest)

    assert install.install_data(data_dir=tmp_path, overwrite=False) == tmp_path
    assert calls == [
        ("gencode-bed", tmp_path, False),
        ("gencode-feature", tmp_path, False),
        ("blacklists", tmp_path, False),
        ("cgi", tmp_path, False),
        ("manifest", tmp_path, None),
    ]


def test_download_gencode_feature_accepts_short_output_dir(monkeypatch, tmp_path):
    called = {}

    def fake_download(species, version, output_dir, **kwargs):
        called["args"] = (species, version, output_dir)
        called["kwargs"] = kwargs
        return {"promoter.up": tmp_path / "promoter.up.bed"}

    monkeypatch.setattr(cli, "download_gencode_feature", fake_download)

    assert cli.main(
        [
            "download-feature",
            "hg38",
            "v31",
            "-o",
            str(tmp_path),
            "-e",
            "3kb",
        ]
    ) == 0
    assert called["args"] == ("hg38", "v31", str(tmp_path))
    assert called["kwargs"]["tes_bp"] == "3kb"


def test_download_gencode_feature_accepts_all_defaults(monkeypatch, tmp_path):
    calls = []

    def fake_download(species, version, output_dir, **kwargs):
        calls.append((species, version, output_dir, kwargs))
        return {"list": tmp_path / species / version / "order.lst"}

    monkeypatch.setattr(cli, "download_gencode_feature", fake_download)

    assert (
        cli.main(["download-feature", "all", "all", "-o", str(tmp_path)])
        == 0
    )
    assert calls
    assert calls[0][0:3] == (
        "hg38",
        "v31",
        str(tmp_path / "hg38" / "v31" / "2kb"),
    )


def test_install_gencode_feature_accepts_feature_options(monkeypatch, tmp_path):
    called = {}

    def fake_install(data_dir=None, overwrite=True, **kwargs):
        called["data_dir"] = data_dir
        called["overwrite"] = overwrite
        called.update(kwargs)
        return tmp_path

    monkeypatch.setattr(cli, "install_gencode_features", fake_install)

    assert (
        cli.main(
            [
                "install-feature",
                "--species",
                "hg38",
                "--version",
                "v31",
                "-d",
                str(tmp_path),
                "-g",
                "genes.gtf",
                "-P",
                "custom",
                "-D",
                "40kb",
                "-n",
            ]
        )
        == 0
    )
    assert called["data_dir"] == str(tmp_path)
    assert called["overwrite"] is False
    assert called["species"] == "hg38"
    assert called["version"] == "v31"
    assert called["gtf_path"] == "genes.gtf"
    assert called["prefix"] == "custom"
    assert called["distal_bp"] == "40kb"


def test_install_gencode_feature_accepts_positional_scope_and_output_dir(
    monkeypatch, tmp_path
):
    called = {}

    def fake_install(data_dir=None, overwrite=True, **kwargs):
        called["data_dir"] = data_dir
        called["overwrite"] = overwrite
        called.update(kwargs)
        return tmp_path

    monkeypatch.setattr(cli, "install_gencode_features", fake_install)

    assert (
        cli.main(
            [
                "install-feature",
                "hg38",
                "all",
                "-d",
                str(tmp_path / "cache"),
                "-o",
                str(tmp_path / "stage"),
            ]
        )
        == 0
    )
    assert called["data_dir"] == str(tmp_path / "cache")
    assert called["overwrite"] is False
    assert called["species"] == "hg38"
    assert called["version"] == "all"
    assert called["output_dir"] == str(tmp_path / "stage")


def test_install_gencode_feature_accepts_all_defaults(monkeypatch, tmp_path):
    called = {}

    def fake_install(data_dir=None, overwrite=True, **kwargs):
        called["data_dir"] = data_dir
        called["overwrite"] = overwrite
        called.update(kwargs)
        return tmp_path

    monkeypatch.setattr(cli, "install_gencode_features", fake_install)

    assert cli.main(["install-feature", "--all", "-d", str(tmp_path)]) == 0
    assert called["data_dir"] == str(tmp_path)
    assert called["overwrite"] is False
    assert called["species"] is None
    assert called["version"] is None


def test_install_gencode_feature_accepts_explicit_overwrite(monkeypatch, tmp_path):
    called = {}

    def fake_install(data_dir=None, overwrite=True, **kwargs):
        called["data_dir"] = data_dir
        called["overwrite"] = overwrite
        called.update(kwargs)
        return tmp_path

    monkeypatch.setattr(cli, "install_gencode_features", fake_install)

    assert (
        cli.main(
            [
                "install-feature",
                "--all",
                "-d",
                str(tmp_path),
                "--overwrite",
            ]
        )
        == 0
    )
    assert called["data_dir"] == str(tmp_path)
    assert called["overwrite"] is True
    assert called["species"] is None
    assert called["version"] is None


def test_install_gencode_feature_requires_scope(capsys):
    assert cli.main(["install-feature"]) == 2
    captured = capsys.readouterr()
    assert "requires species and version, or --all" in captured.err


def test_download_chromhmm_accepts_short_options(monkeypatch, tmp_path):
    called = {}

    def fake_download(
        data_dir=None,
        model=18,
        genome="hg19",
        ids=None,
        tissue=None,
        cellline=None,
        overwrite=True,
        progress=None,
    ):
        called["data_dir"] = data_dir
        called["model"] = model
        called["genome"] = genome
        called["ids"] = ids
        called["tissue"] = tissue
        called["cellline"] = cellline
        called["overwrite"] = overwrite
        called["progress"] = progress
        return {"E063": tmp_path / "E063.bed.gz"}

    monkeypatch.setattr(cli, "download_chromhmm", fake_download)

    assert (
        cli.main(
            [
                "download-chromhmm",
                "-o",
                str(tmp_path),
                "-m",
                "15",
                "-G",
                "hg38",
                "-i",
                "E063",
                "-t",
                "blood",
                "-c",
                "GM12878",
                "-n",
            ]
        )
        == 0
    )
    assert called == {
        "data_dir": str(tmp_path),
        "model": "15",
        "genome": "hg38",
        "ids": ["E063"],
        "tissue": ["blood"],
        "cellline": ["GM12878"],
        "overwrite": False,
        "progress": cli._stderr_progress,
    }


def test_install_chromhmm_uses_data_dir_option(monkeypatch, tmp_path):
    called = {}

    def fake_download(
        data_dir=None,
        model=18,
        genome="hg19",
        ids=None,
        tissue=None,
        cellline=None,
        overwrite=True,
        progress=None,
    ):
        called["data_dir"] = data_dir
        called["model"] = model
        called["genome"] = genome
        called["ids"] = ids
        called["tissue"] = tissue
        called["cellline"] = cellline
        called["overwrite"] = overwrite
        called["progress"] = progress
        return {"E063": tmp_path / "E063.bed.gz"}

    monkeypatch.setattr(cli, "download_chromhmm", fake_download)

    assert (
        cli.main(["install-chromhmm", "-d", str(tmp_path), "-i", "E063", "-n"])
        == 0
    )
    assert called["data_dir"] == str(tmp_path)
    assert called["ids"] == ["E063"]
    assert called["overwrite"] is False


def test_download_chromhmm_by_id_builds_dense_bed_urls(monkeypatch, tmp_path):
    downloaded = {}

    def fake_download(url, destination, progress=None):
        downloaded[destination.name] = url
        destination.write_text(url, encoding="utf-8")

    def fake_fetch(url):
        assert url == chromhmm.CHROMHMM_METADATA_URL
        return _chromhmm_data_js()

    monkeypatch.setattr(
        chromhmm,
        "chromhmm_available_ids",
        lambda model, genome="hg19": ("E003", "E063"),
    )
    monkeypatch.setattr(chromhmm, "_fetch_text", fake_fetch)
    monkeypatch.setattr(chromhmm, "_download_file", fake_download)

    outputs = db.download_chromhmm(data_dir=tmp_path, model=18, ids="E003,E063")

    assert outputs["E003"] == tmp_path / "chromhmm" / "hg19" / "18state" / (
        "E003_18_core_K27ac_dense.bed.gz"
    )
    assert downloaded["E003_18_core_K27ac_dense.bed.gz"] == db.chromhmm_url(
        "E003", 18
    )
    assert downloaded["E063_18_core_K27ac_dense.bed.gz"] == db.chromhmm_url(
        "E063", 18
    )
    log_text = (tmp_path / "download_urls.log").read_text(encoding="utf-8")
    assert db.chromhmm_url("E003", 18) in log_text
    assert db.chromhmm_url("E063", 18) in log_text
    assert chromhmm.CHROMHMM_METADATA_URL in log_text
    metadata_tsv = tmp_path / "chromhmm" / "metadata.tsv"
    assert metadata_tsv.read_text(encoding="utf-8") == (
        "path\tversion\tstate_number\teid\tgroup\tmnemonic\tname\t"
        "name_edacc_r9\tcategory\n"
        "{}\thg19\t18\tE003\tBlood_T_cell\tBLD_CD19\t"
        "GM12878_Lymphoblastoid_Cells\t"
        "GM12878_Lymphoblastoid_Cell_Line\tClass_1\n"
        "{}\thg19\t18\tE063\tBrain_Fetal\tBRN_FETAL\tFetal_Brain_Male\t"
        "Fetal_Brain_Male\tClass_2\n"
    ).format(
        tmp_path / "chromhmm" / "hg19" / "18state" / (
            "E003_18_core_K27ac_dense.bed.gz"
        ),
        tmp_path / "chromhmm" / "hg19" / "18state" / (
            "E063_18_core_K27ac_dense.bed.gz"
        ),
    )


def test_download_chromhmm_skips_unavailable_requested_ids(monkeypatch, tmp_path):
    downloaded = {}

    def fake_download(url, destination, progress=None):
        downloaded[destination.name] = url
        destination.write_text(url, encoding="utf-8")

    def fake_fetch(url):
        assert url == chromhmm.CHROMHMM_METADATA_URL
        return _chromhmm_data_js()

    monkeypatch.setattr(
        chromhmm,
        "chromhmm_available_ids",
        lambda model, genome="hg19": ("E063",),
    )
    monkeypatch.setattr(chromhmm, "_fetch_text", fake_fetch)
    monkeypatch.setattr(chromhmm, "_download_file", fake_download)

    outputs = db.download_chromhmm(data_dir=tmp_path, model=18, ids="E001,E063")

    assert tuple(outputs) == ("E063",)
    assert tuple(downloaded) == ("E063_18_core_K27ac_dense.bed.gz",)


def test_download_chromhmm_hg38_uses_lifted_over_url(monkeypatch, tmp_path):
    downloaded = {}

    def fake_download(url, destination, progress=None):
        downloaded[destination.name] = url
        destination.write_text(url, encoding="utf-8")

    def fake_fetch(url):
        assert url == chromhmm.CHROMHMM_METADATA_URL
        return _chromhmm_data_js()

    monkeypatch.setattr(
        chromhmm,
        "chromhmm_available_ids",
        lambda model, genome="hg19": ("E063",),
    )
    monkeypatch.setattr(chromhmm, "_fetch_text", fake_fetch)
    monkeypatch.setattr(chromhmm, "_download_file", fake_download)

    outputs = db.download_chromhmm(
        data_dir=tmp_path,
        model=18,
        genome="hg38",
        ids="E063",
    )

    assert outputs["E063"] == tmp_path / "chromhmm" / "hg38" / "18state" / (
        "E063_18_core_K27ac_dense.bed.gz"
    )
    assert downloaded["E063_18_core_K27ac_dense.bed.gz"] == (
        "https://egg2.wustl.edu/roadmap/data/byFileType/"
        "chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/"
        "bed_hg38_lifted_over/E063_18_core_K27ac_dense.bed.gz"
    )


def test_download_chromhmm_fuzzy_cellline_selection(monkeypatch, tmp_path):
    data_js = """
var data_epg=
[
  {
    "eid":"E001",
    "group":"ESC",
    "mnemonic":"ESC.I3",
    "name":"ES-I3 Cells",
    "name_EDACC_R9":"ES-I3_Cell_Line",
    "class":"Class 4",
    "info":""
  },
  {
    "eid":"E063",
    "group":"Blood & T-cell",
    "mnemonic":"BLD.CD19",
    "name":"GM12878 Lymphoblastoid Cells",
    "name_EDACC_R9":"GM12878_Lymphoblastoid_Cell_Line",
    "class":"Class 1",
    "info":""
  }
];
"""

    def fake_fetch(url):
        if url == chromhmm.CHROMHMM_METADATA_URL:
            return data_js
        return '<a href="E063_18_core_K27ac_dense.bed.gz">dense</a>'

    def fake_download(url, destination):
        destination.write_text(url, encoding="utf-8")

    monkeypatch.setattr(chromhmm, "_fetch_text", fake_fetch)
    monkeypatch.setattr(chromhmm, "_download_file", fake_download)

    outputs = db.download_chromhmm(
        data_dir=tmp_path,
        model=18,
        cellline="gm12878",
    )

    assert tuple(outputs) == ("E063",)
    assert outputs["E063"].name == "E063_18_core_K27ac_dense.bed.gz"


def test_download_segway_accepts_short_options(monkeypatch, tmp_path):
    called = {}

    def fake_download(
        data_dir=None,
        genome="hg19",
        names=None,
        tissue=None,
        cellline=None,
        all_celltypes=False,
        include_encyclopedia=False,
        include_caas=False,
        include_label_info=False,
        overwrite=True,
        progress=None,
    ):
        called["data_dir"] = data_dir
        called["genome"] = genome
        called["names"] = names
        called["tissue"] = tissue
        called["cellline"] = cellline
        called["all_celltypes"] = all_celltypes
        called["include_encyclopedia"] = include_encyclopedia
        called["include_caas"] = include_caas
        called["include_label_info"] = include_label_info
        called["overwrite"] = overwrite
        called["progress"] = progress
        return {"GM12878": tmp_path / "GM12878.bed.gz"}

    monkeypatch.setattr(cli, "download_segway", fake_download)

    assert (
        cli.main(
            [
                "download-segway",
                "-o",
                str(tmp_path),
                "-G",
                "hg19",
                "-i",
                "GM12878",
                "-t",
                "brain",
                "-c",
                "H1",
                "--all-celltypes",
                "--include-encyclopedia",
                "--include-caas",
                "--include-label-info",
            ]
        )
        == 0
    )
    assert called == {
        "data_dir": str(tmp_path),
        "genome": "hg19",
        "names": ["GM12878"],
        "tissue": ["brain"],
        "cellline": ["H1"],
        "all_celltypes": True,
        "include_encyclopedia": True,
        "include_caas": True,
        "include_label_info": True,
        "overwrite": False,
        "progress": cli._stderr_progress,
    }


def test_install_segway_uses_data_dir_option(monkeypatch, tmp_path):
    called = {}

    def fake_install(
        data_dir=None,
        genome="hg19",
        names=None,
        tissue=None,
        cellline=None,
        all_celltypes=False,
        include_encyclopedia=False,
        include_caas=False,
        include_label_info=False,
        overwrite=True,
        progress=None,
    ):
        called["data_dir"] = data_dir
        called["names"] = names
        called["overwrite"] = overwrite
        return {"GM12878": tmp_path / "GM12878.bed.gz"}

    monkeypatch.setattr(cli, "install_segway", fake_install)

    assert cli.main(["install-segway", "-d", str(tmp_path), "-i", "GM12878", "-n"]) == 0
    assert called["data_dir"] == str(tmp_path)
    assert called["names"] == ["GM12878"]
    assert called["overwrite"] is False


def test_download_segway_discovers_encode_accessions(monkeypatch, tmp_path):
    downloaded = {}
    messages = []
    source_url = segway.SEGWAY_URL.replace("www.", "test.")

    def fake_fetch_json(url):
        assert url == (
            source_url
            + "?format=json&frame=embedded"
        )
        return {
            "@id": "/publications/94941f71-80c8-43d2-809b-25161efc3be0/",
            "files": [
                _segway_encode_file(
                    "ENCFF314OLC",
                    "/files/ENCFF314OLC/@@download/ENCFF314OLC.bed.gz",
                    "segway_encyclopedia.bed.gz",
                    output_type="Segway encyclopedia",
                ),
                _segway_encode_file(
                    "ENCFF222ABC",
                    "/files/ENCFF222ABC/@@download/ENCFF222ABC.bed.gz",
                    (
                        "gs://encode/segway/"
                        "PERIPHERAL_BLOOD_MONONUCLEAR_PRIMARY_CELLS/ann.bed"
                    ),
                ),
            ],
        }

    def fake_download(url, destination, progress=None):
        downloaded[destination.name] = url
        destination.write_text(url, encoding="utf-8")

    monkeypatch.setattr(segway, "_fetch_encode_json", fake_fetch_json)
    monkeypatch.setattr(segway, "_download_file", fake_download)

    outputs = db.download_segway(
        data_dir=tmp_path,
        source_url=source_url,
        fallback_url=source_url,
        names="PERIPHERAL_BLOOD_MONONUCLEAR_PRIMARY_CELLS",
        include_encyclopedia=True,
        progress=messages.append,
    )

    assert outputs["PERIPHERAL_BLOOD_MONONUCLEAR_PRIMARY_CELLS"] == (
        tmp_path / "segway" / "hg19" / "PERIPHERAL_BLOOD_MONONUCLEAR_PRIMARY_CELLS.bed.gz"
    )
    assert outputs["segway_encyclopedia"] == (
        tmp_path / "segway" / "hg19" / "segway_encyclopedia.bed.gz"
    )
    assert downloaded == {
        "PERIPHERAL_BLOOD_MONONUCLEAR_PRIMARY_CELLS.bed.gz": (
            "https://www.encodeproject.org/files/ENCFF222ABC/@@download/"
            "ENCFF222ABC.bed.gz"
        ),
        "segway_encyclopedia.bed.gz": (
            "https://www.encodeproject.org/files/ENCFF314OLC/@@download/"
            "ENCFF314OLC.bed.gz"
        ),
    }
    log_text = (tmp_path / "download_urls.log").read_text(encoding="utf-8")
    assert "ENCFF314OLC" in log_text
    assert "ENCFF222ABC" in log_text
    assert not any("fallback" in message for message in messages)


def test_download_segway_searches_encode_related_accessions(monkeypatch, tmp_path):
    downloaded = {}
    seen_urls = []
    source_url = segway.SEGWAY_URL.replace("www.", "test.")

    def fake_fetch_json(url):
        seen_urls.append(url)
        if url == source_url + "?format=json&frame=embedded":
            return {
                "@id": "/publications/94941f71-80c8-43d2-809b-25161efc3be0/",
                "files": [],
            }
        if "type=File" in url and "references=" in url:
            return {
                "@graph": [
                    _segway_encode_file(
                        "ENCFF555AAA",
                        "/files/ENCFF555AAA/@@download/ENCFF555AAA.bed.gz",
                        "gs://encode/segway/H1-HESC/ann.bed",
                    )
                ]
            }
        if "type=Annotation" in url and "references=" in url:
            return {"@graph": []}
        raise AssertionError(url)

    def fake_download(url, destination, progress=None):
        downloaded[destination.name] = url
        destination.write_text(url, encoding="utf-8")

    monkeypatch.setattr(segway, "_fetch_encode_json", fake_fetch_json)
    monkeypatch.setattr(segway, "_download_file", fake_download)

    outputs = db.download_segway(
        data_dir=tmp_path,
        source_url=source_url,
        fallback_url=source_url,
        names="H1-HESC",
    )

    assert outputs["H1-HESC"] == tmp_path / "segway" / "hg19" / "H1-HESC.bed.gz"
    assert downloaded == {
        "H1-HESC.bed.gz": (
            "https://www.encodeproject.org/files/ENCFF555AAA/@@download/"
            "ENCFF555AAA.bed.gz"
        )
    }
    assert any("type=File" in url for url in seen_urls)


def test_download_segway_uses_bundled_manifest_for_blood(monkeypatch, tmp_path):
    downloaded = {}
    messages = []

    def fail_fetch(*args, **kwargs):
        raise AssertionError("built-in Segway sources should use bundled metadata")

    def fake_download(url, destination, progress=None):
        downloaded[destination.name] = url
        destination.write_text(url, encoding="utf-8")

    monkeypatch.setattr(segway, "_fetch_encode_json", fail_fetch)
    monkeypatch.setattr(segway, "_fetch_text", fail_fetch)
    monkeypatch.setattr(segway, "_download_file", fake_download)
    monkeypatch.setattr(segway, "_BUNDLED_SEGWAY_FILES", None)

    outputs = db.download_segway(
        data_dir=tmp_path,
        tissue="blood",
        progress=messages.append,
    )

    assert outputs == {
        "PERIPHERAL_BLOOD_MONONUCLEAR_PRIMARY_CELLS": (
            tmp_path
            / "segway"
            / "hg19"
            / "PERIPHERAL_BLOOD_MONONUCLEAR_PRIMARY_CELLS.bed.gz"
        )
    }
    assert downloaded == {
        "PERIPHERAL_BLOOD_MONONUCLEAR_PRIMARY_CELLS.bed.gz": (
            "https://www.encodeproject.org/files/ENCFF659DOR/@@download/ENCFF659DOR.bed.gz"
        )
    }
    db.download_segway(data_dir=tmp_path, tissue="blood", progress=messages.append)
    metadata_lines = (tmp_path / "segway" / "metadata.tsv").read_text(
        encoding="utf-8"
    ).splitlines()
    assert (
        sum(
            "PERIPHERAL_BLOOD_MONONUCLEAR_PRIMARY_CELLS" in line
            for line in metadata_lines
        )
        == 1
    )
    assert not any("primary source unavailable" in message for message in messages)
    assert not any("primary source had no matching" in message for message in messages)
    manifest_text = db.segway_download_manifest_path(tmp_path).read_text(
        encoding="utf-8"
    )
    assert "ENCFF338HEJ" in manifest_text
    assert "PERIPHERAL_BLOOD_MONONUCLEAR_PRIMARY_CELLS" in manifest_text


def test_download_segway_uses_bundled_encode_manifest(monkeypatch, tmp_path):
    downloaded = {}
    messages = []

    def fail_fetch(*args, **kwargs):
        raise AssertionError("built-in Segway sources should use bundled metadata")

    def fake_download(url, destination, progress=None):
        downloaded[destination.name] = url
        destination.write_text(url, encoding="utf-8")

    monkeypatch.setattr(segway, "_fetch_encode_json", fail_fetch)
    monkeypatch.setattr(segway, "_fetch_text", fail_fetch)
    monkeypatch.setattr(segway, "_download_file", fake_download)
    monkeypatch.setattr(segway, "_BUNDLED_SEGWAY_FILES", None)

    outputs = db.download_segway(
        data_dir=tmp_path,
        names="COLONIC_MUCOSA",
        progress=messages.append,
    )

    encode_url = (
        "https://www.encodeproject.org/files/ENCFF725HYN/@@download/"
        "ENCFF725HYN.bed.gz"
    )
    assert outputs == {
        "COLONIC_MUCOSA": tmp_path
        / "segway"
        / "hg19"
        / "COLONIC_MUCOSA.bed.gz"
    }
    assert downloaded == {"COLONIC_MUCOSA.bed.gz": encode_url}
    assert not any("fallback" in message for message in messages)
    assert encode_url in (tmp_path / "segway" / "metadata.tsv").read_text(
        encoding="utf-8"
    )
    assert "ENCFF338HEJ" in db.segway_download_manifest_path(tmp_path).read_text(
        encoding="utf-8"
    )


def test_segway_bundled_manifest_combines_encode_and_washington(monkeypatch):
    def fail_fetch(*args, **kwargs):
        raise AssertionError("built-in Segway sources should use bundled metadata")

    monkeypatch.setattr(segway, "_fetch_encode_json", fail_fetch)
    monkeypatch.setattr(segway, "_fetch_text", fail_fetch)
    monkeypatch.setattr(segway, "_BUNDLED_SEGWAY_FILES", None)

    names = db.segway_available_names()

    assert "COLONIC_MUCOSA" in names
    assert "PERIPHERAL_BLOOD_MONONUCLEAR_PRIMARY_CELLS" in names
    assert db.segway_url("COLONIC_MUCOSA") == (
        "https://www.encodeproject.org/files/ENCFF725HYN/@@download/"
        "ENCFF725HYN.bed.gz"
    )
    assert db.segway_url("PERIPHERAL_BLOOD_MONONUCLEAR_PRIMARY_CELLS") == (
        "https://www.encodeproject.org/files/ENCFF659DOR/@@download/"
        "ENCFF659DOR.bed.gz"
    )


def test_download_segway_by_name_builds_urls(monkeypatch, tmp_path):
    downloaded = {}
    source_url = "https://example.org/segway/"
    interpreted_url = source_url + "interpreted/"

    def fake_fetch(url):
        if url == source_url:
            return _segway_index_html()
        if url == interpreted_url:
            return _segway_interpreted_html()
        raise AssertionError(url)

    def fake_download(url, destination, progress=None):
        downloaded[destination.name] = url
        destination.write_text(url, encoding="utf-8")

    monkeypatch.setattr(segway, "_fetch_text", fake_fetch)
    monkeypatch.setattr(segway, "_download_file", fake_download)

    outputs = db.download_segway(
        data_dir=tmp_path,
        names="GM12878,H1-HESC",
        source_url=source_url,
        fallback_url=source_url,
    )

    assert outputs["GM12878"] == (
        tmp_path / "segway" / "hg19" / "GM12878.bed.gz"
    )
    assert outputs["H1-HESC"] == (
        tmp_path / "segway" / "hg19" / "H1-HESC.bed.gz"
    )
    assert downloaded == {
        "GM12878.bed.gz": (
            interpreted_url + "GM12878.bed.gz"
        ),
        "H1-HESC.bed.gz": (
            interpreted_url + "H1-HESC.bed.gz"
        ),
    }
    metadata_text = (tmp_path / "segway" / "metadata.tsv").read_text(
        encoding="utf-8"
    )
    assert "encodeID\tbiosample_term_name" in metadata_text
    assert "GM12878\tcelltype" in metadata_text
    assert "H1-HESC\tcelltype" in metadata_text
    log_text = (tmp_path / "download_urls.log").read_text(encoding="utf-8")
    assert "interpreted/GM12878.bed.gz" in log_text
    assert source_url in log_text


def test_download_segway_skips_existing_files_by_default(monkeypatch, tmp_path):
    destination = tmp_path / "segway" / "hg19" / "GM12878.bed.gz"
    destination.parent.mkdir(parents=True)
    destination.write_text("cached\n", encoding="utf-8")
    source_url = "https://example.org/segway/"
    interpreted_url = source_url + "interpreted/"

    def fake_fetch(url):
        if url == source_url:
            return _segway_index_html()
        if url == interpreted_url:
            return _segway_interpreted_html()
        raise AssertionError(url)

    def fake_download(url, destination):
        raise AssertionError("existing Segway files should be skipped")

    monkeypatch.setattr(segway, "_fetch_text", fake_fetch)
    monkeypatch.setattr(segway, "_download_file", fake_download)

    outputs = db.download_segway(
        data_dir=tmp_path,
        names="GM12878",
        source_url=source_url,
        fallback_url=source_url,
    )

    assert outputs["GM12878"] == destination
    assert destination.read_text(encoding="utf-8") == "cached\n"


def test_download_segway_defaults_to_encyclopedia(monkeypatch, tmp_path):
    downloaded = {}
    source_url = "https://example.org/segway/"
    interpreted_url = source_url + "interpreted/"

    def fake_fetch(url):
        if url == source_url:
            return _segway_index_html()
        if url == interpreted_url:
            return _segway_interpreted_html()
        raise AssertionError(url)

    def fake_download(url, destination):
        downloaded[destination.name] = url
        destination.write_text(url, encoding="utf-8")

    monkeypatch.setattr(segway, "_fetch_text", fake_fetch)
    monkeypatch.setattr(segway, "_download_file", fake_download)

    outputs = db.download_segway(
        data_dir=tmp_path,
        source_url=source_url,
        fallback_url=source_url,
    )

    assert tuple(outputs) == ("segway_encyclopedia",)
    assert outputs["segway_encyclopedia"] == (
        tmp_path / "segway" / "hg19" / "segway_encyclopedia.bed.gz"
    )
    assert downloaded == {
        "segway_encyclopedia.bed.gz": (
            source_url + "segway_encyclopedia.bed.gz"
        )
    }


def test_download_segway_fuzzy_cellline_selection(monkeypatch, tmp_path):
    source_url = "https://example.org/segway/"
    interpreted_url = source_url + "interpreted/"

    def fake_fetch(url):
        if url == source_url:
            return _segway_index_html()
        if url == interpreted_url:
            return _segway_interpreted_html()
        raise AssertionError(url)

    def fake_download(url, destination):
        destination.write_text(url, encoding="utf-8")

    monkeypatch.setattr(segway, "_fetch_text", fake_fetch)
    monkeypatch.setattr(segway, "_download_file", fake_download)

    outputs = db.download_segway(
        data_dir=tmp_path,
        cellline="h1 hesc",
        source_url=source_url,
        fallback_url=source_url,
    )

    assert tuple(outputs) == ("H1-HESC",)
    assert outputs["H1-HESC"].name == "H1-HESC.bed.gz"


def test_download_segway_uses_fallback_source(monkeypatch, tmp_path):
    downloaded = {}
    source_url = segway.SEGWAY_URL.replace("www.", "test.")
    fallback_url = "https://fallback.example.org/segway/"
    messages = []

    def fake_fetch_json(url):
        raise ValueError("human verification/captcha page returned")

    def fake_fetch(url):
        if url == fallback_url:
            return _segway_index_html()
        if url == fallback_url + "interpreted/":
            return _segway_interpreted_html()
        raise AssertionError(url)

    def fake_download(url, destination, progress=None):
        downloaded[destination.name] = url
        destination.write_text(url, encoding="utf-8")

    monkeypatch.setattr(segway, "_fetch_encode_json", fake_fetch_json)
    monkeypatch.setattr(segway, "_fetch_text", fake_fetch)
    monkeypatch.setattr(segway, "_download_file", fake_download)

    outputs = db.download_segway(
        data_dir=tmp_path,
        source_url=source_url,
        fallback_url=fallback_url,
        progress=messages.append,
    )

    assert tuple(outputs) == ("segway_encyclopedia",)
    assert downloaded == {
        "segway_encyclopedia.bed.gz": (
            fallback_url + "segway_encyclopedia.bed.gz"
        )
    }
    assert fallback_url in (tmp_path / "download_urls.log").read_text(
        encoding="utf-8"
    )
    assert any("primary source unavailable" in message for message in messages)
    assert any("using fallback source" in message for message in messages)


def test_download_segway_non_hg19_requires_liftover_confirmation(capsys):
    assert cli.main(["download-segway", "-G", "hg38"]) == 2
    captured = capsys.readouterr()
    assert "--yes-liftover" in captured.err


def test_download_segway_yes_liftover_writes_script(monkeypatch, tmp_path, capsys):
    called = {}

    def fake_download(
        data_dir=None,
        genome="hg19",
        names=None,
        tissue=None,
        cellline=None,
        all_celltypes=False,
        include_encyclopedia=False,
        include_caas=False,
        include_label_info=False,
        overwrite=True,
        progress=None,
    ):
        called["data_dir"] = data_dir
        called["genome"] = genome
        called["names"] = names
        return {"GM12878": tmp_path / "GM12878.bed.gz"}

    def fake_write(
        data_dir=None,
        target_genome="hg38",
        script_path=None,
        install_root=None,
    ):
        called["liftover_data_dir"] = data_dir
        called["target_genome"] = target_genome
        called["install_root"] = install_root
        return tmp_path / "segway" / "liftover_hg19_to_hg38.sh"

    monkeypatch.setattr(cli, "download_segway", fake_download)
    monkeypatch.setattr(cli, "write_segway_liftover_script", fake_write)

    assert (
        cli.main(
            [
                "download-segway",
                "-o",
                str(tmp_path),
                "-G",
                "hg38",
                "--yes-liftover",
                "-i",
                "GM12878",
            ]
        )
        == 0
    )
    captured = capsys.readouterr()
    assert called["data_dir"] == str(tmp_path)
    assert called["genome"] == "hg19"
    assert called["target_genome"] == "hg38"
    assert called["install_root"] is None
    assert "liftover_script\t" in captured.out


def test_install_segway_yes_liftover_sets_install_root(monkeypatch, tmp_path):
    called = {}

    def fake_install(
        data_dir=None,
        genome="hg19",
        names=None,
        tissue=None,
        cellline=None,
        all_celltypes=False,
        include_encyclopedia=False,
        include_caas=False,
        include_label_info=False,
        overwrite=False,
        progress=None,
    ):
        called["data_dir"] = data_dir
        called["genome"] = genome
        return {"GM12878": tmp_path / "GM12878.bed.gz"}

    def fake_write(
        data_dir=None,
        target_genome="hg38",
        script_path=None,
        install_root=None,
    ):
        called["liftover_data_dir"] = data_dir
        called["target_genome"] = target_genome
        called["install_root"] = install_root
        return tmp_path / "segway" / "liftover_hg19_to_hg38.sh"

    monkeypatch.setattr(cli, "install_segway", fake_install)
    monkeypatch.setattr(cli, "write_segway_liftover_script", fake_write)

    assert (
        cli.main(
            [
                "install-segway",
                "-d",
                str(tmp_path),
                "-G",
                "hg38",
                "--yes-liftover",
                "-i",
                "GM12878",
            ]
        )
        == 0
    )

    assert called["data_dir"] == str(tmp_path)
    assert called["genome"] == "hg19"
    assert called["liftover_data_dir"] == str(tmp_path)
    assert called["target_genome"] == "hg38"
    assert called["install_root"] == tmp_path / "segway"


def test_segway_liftover_script_uses_crossmap_and_ucsc_chain(tmp_path):
    script = db.segway_liftover_script(
        tmp_path / "segway" / "hg19",
        tmp_path / "segway" / "hg38",
        "hg38",
        install_root=tmp_path / "install" / "segway",
    )

    assert "env list" not in script
    assert "grep -Fxq" not in script
    assert "activate_existing_env" in script
    assert "Activated existing conda environment" in script
    assert "Activated existing micromamba environment" in script
    assert "Activated existing mamba environment" in script
    assert "micromamba activate" in script
    assert "conda activate" in script
    assert "mamba activate" in script
    assert "conda run" not in script
    assert "install -y -n" not in script
    assert "remove the env so this script can recreate it" in script
    assert "create_and_activate_env" in script
    assert 'create -y -n "${ENV_NAME}"' in script
    assert "crossmap" in script
    assert "cp -p" in script
    assert "INSTALL_ROOT=" in script
    assert "hg38lift.bed" in script
    assert '${OUTPUT_DIR}"/interpreted' not in script
    assert '${INSTALL_DIR}"/interpreted' not in script
    assert "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/" in script
    assert "hg19ToHg38.over.chain.gz" in script
    assert "CrossMap" in script


def test_write_gencode_region_unions_from_gtf(tmp_path):
    gtf = _write_mini_gencode_region_gtf(tmp_path)

    outputs = db.write_gencode_region_unions(gtf, tmp_path, tes_bp=2)

    assert outputs["exon"] == tmp_path / "exon.bed"
    assert outputs["exon"].read_text(encoding="utf-8") == (
        "chr1\t100\t150\n"
        "chr1\t180\t220\n"
        "chr1\t400\t450\n"
        "chr1\t480\t500\n"
    )
    assert outputs["intron"].read_text(encoding="utf-8") == (
        "chr1\t120\t200\n"
        "chr1\t450\t480\n"
    )
    assert outputs["tes"].read_text(encoding="utf-8") == (
        "chr1\t217\t222\n"
        "chr1\t398\t403\n"
    )
    assert outputs["intergenic"].read_text(encoding="utf-8") == (
        "chr1\t0\t100\n"
        "chr1\t220\t400\n"
        "chr1\t500\t600\n"
    )


def test_download_gencode_regions_writes_gtf_derived_beds(tmp_path):
    gtf = _write_mini_gencode_region_gtf(tmp_path)

    outputs = db.download_gencode_tss_flank_region_unions(
        "hg38",
        "v31",
        tmp_path,
        gtf_path=gtf,
        promoter_bp=20,
        distal_bp=50,
        tes_bp=2,
    )

    assert outputs["gene_bed"] == tmp_path / "gencode.v31.hg38.gene.bed.withtype"
    assert outputs["exon"] == tmp_path / "20bp.exon.bed"
    assert outputs["intron"] == tmp_path / "20bp.intron.bed"
    assert outputs["tes"] == tmp_path / "20bp.tes.bed"
    assert outputs["intergenic"] == tmp_path / "20bp.intergenic.bed"
    assert outputs["promoter"] == tmp_path / "20bp.promoter.bed"
    assert outputs["list"] == tmp_path / "order.lst"
    assert outputs["list"].read_text(encoding="utf-8") == (
        "20bp.promoter.up.bed\n"
        "20bp.promoter.down.bed\n"
        "20bp.exon.bed\n"
        "20bp.intron.bed\n"
        "20bp.tes.bed\n"
        "20bp.dis5.bed\n"
        "20bp.dis3.bed\n"
        "20bp.intergenic.bed\n"
    )
    assert outputs["tes"].read_text(encoding="utf-8") == (
        "chr1\t217\t222\n"
        "chr1\t398\t403\n"
    )
    assert outputs["promoter"].read_text(encoding="utf-8") == (
        "chr1\t80\t121\n"
        "chr1\t479\t520\n"
    )
    assert outputs["dis3"].read_text(encoding="utf-8") == (
        "chr1\t240\t270\n"
        "chr1\t350\t380\n"
    )
    assert outputs["intergenic"].read_text(encoding="utf-8") == (
        "chr1\t0\t50\n"
        "chr1\t222\t240\n"
        "chr1\t270\t350\n"
        "chr1\t380\t398\n"
        "chr1\t550\t600\n"
    )


def test_install_gencode_feature_set_uses_feature_prefix_layout(tmp_path):
    gtf = _write_mini_gencode_region_gtf(tmp_path)

    feature_dir = install.install_gencode_feature_set(
        "hg38",
        "v31",
        data_dir=tmp_path,
        gtf_path=gtf,
        promoter_bp="2kb",
        distal_bp="50kb",
        tes_bp=2,
    )

    version_dir = tmp_path / "feature" / "hg38" / "v31"
    assert feature_dir == version_dir / "2kb"
    assert (feature_dir / "2kb.promoter.up.bed").exists()
    assert (feature_dir / "2kb.exon.bed").exists()
    assert (feature_dir / "2kb.intron.bed").exists()
    assert (feature_dir / "2kb.tes.bed").exists()
    assert (feature_dir / "2kb.intergenic.bed").exists()
    assert (feature_dir / "2kb.promoter.bed").exists()
    assert not (feature_dir / "gencode.v31.hg38.gene.bed.withtype").exists()
    assert (tmp_path / "bed" / "hg38" / "v31" / "all.gene.bed").exists()
    assert (feature_dir / "order.lst").read_text(encoding="utf-8") == (
        "2kb.promoter.up.bed\n"
        "2kb.promoter.down.bed\n"
        "2kb.exon.bed\n"
        "2kb.intron.bed\n"
        "2kb.tes.bed\n"
        "2kb.dis5.bed\n"
        "2kb.dis3.bed\n"
        "2kb.intergenic.bed\n"
    )
    default_dir = tmp_path / "feature" / "hg38" / "def"
    if default_dir.is_symlink():
        assert os.readlink(str(default_dir)) == "v31/2kb"
    assert (default_dir / "order.lst").exists()


def test_install_gencode_feature_set_reuses_preprocessed_output(monkeypatch, tmp_path):
    source_dir = tmp_path / "stage"
    label = "2kb"
    source_dir.mkdir()
    order_entries = [
        "{}.{}.bed".format(label, region_type)
        for region_type in db.GENCODE_FEATURE_LIST_ORDER
    ]
    (source_dir / "order.lst").write_text(
        "{}\n".format("\n".join(order_entries)),
        encoding="utf-8",
    )
    for name in order_entries:
        (source_dir / name).write_text("{}\n".format(name), encoding="utf-8")
    (source_dir / "{}.promoter.bed".format(label)).write_text(
        "combined promoter\n",
        encoding="utf-8",
    )

    def fail_download(*args, **kwargs):
        raise AssertionError("download_gencode_feature should not be called")

    monkeypatch.setattr(install, "download_gencode_feature", fail_download)

    feature_dir = install.install_gencode_feature_set(
        "hg38",
        "v31",
        data_dir=tmp_path / "cache",
        source_dir=source_dir,
        promoter_bp="2kb",
    )

    assert feature_dir == tmp_path / "cache" / "feature" / "hg38" / "v31" / "2kb"
    assert (feature_dir / "order.lst").exists()
    assert (feature_dir / "2kb.exon.bed").read_text(encoding="utf-8") == (
        "2kb.exon.bed\n"
    )
    assert (feature_dir / "2kb.promoter.bed").read_text(encoding="utf-8") == (
        "combined promoter\n"
    )
    default_dir = tmp_path / "cache" / "feature" / "hg38" / "def"
    if default_dir.is_symlink():
        assert os.readlink(str(default_dir)) == "v31/2kb"


def test_install_gencode_feature_set_skips_existing_target(monkeypatch, tmp_path):
    label = "2kb"
    feature_dir = tmp_path / "feature" / "hg38" / "v31" / label
    feature_dir.mkdir(parents=True)
    order_entries = [
        "{}.{}.bed".format(label, region_type)
        for region_type in db.GENCODE_FEATURE_LIST_ORDER
    ]
    (feature_dir / "order.lst").write_text(
        "{}\n".format("\n".join(order_entries)),
        encoding="utf-8",
    )
    for name in order_entries:
        (feature_dir / name).write_text("{}\n".format(name), encoding="utf-8")

    def fail_download(*args, **kwargs):
        raise AssertionError("download_gencode_feature should not be called")

    monkeypatch.setattr(install, "download_gencode_feature", fail_download)

    result = install.install_gencode_feature_set(
        "hg38",
        "v31",
        data_dir=tmp_path,
        promoter_bp="2kb",
        skip_existing=True,
    )

    assert result == feature_dir
    default_dir = tmp_path / "feature" / "hg38" / "def"
    if default_dir.is_symlink():
        assert os.readlink(str(default_dir)) == "v31/2kb"


def test_convert_gencode_gtf_to_bed_withtype(tmp_path):
    gtf = tmp_path / "mini.gtf"
    gtf.write_text(
        "\n".join(
            [
                'chr1\tGENCODE\ttranscript\t101\t220\t.\t+\t.\tgene_id "ENSG1.2"; transcript_id "ENST1.3"; gene_name "GENE1"; gene_type "protein_coding";',
                'chr1\tGENCODE\texon\t101\t120\t.\t+\t.\tgene_id "ENSG1.2"; transcript_id "ENST1.3"; gene_name "GENE1"; gene_type "protein_coding";',
                'chr1\tGENCODE\texon\t201\t220\t.\t+\t.\tgene_id "ENSG1.2"; transcript_id "ENST1.3"; gene_name "GENE1"; gene_type "protein_coding";',
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    output = tmp_path / "genes.bed.withtype"

    db.convert_gencode_gtf_to_bed(gtf, output)

    assert output.read_text(encoding="utf-8") == (
        "chr1\t100\t220\tGENE1\t40\t+\tENSG1.2\tENST1.3\t"
        "protein_coding\n"
    )


def test_write_deduplong_can_group_by_symbol_or_ensid(tmp_path):
    gene_bed = tmp_path / "genes.bed"
    gene_bed.write_text(
        "\n".join(
            [
                "chr1\t100\t200\tSameSymbol\t100\t+\tENSG1.1\tTX1.1",
                "chr1\t100\t250\tSameSymbol\t150\t+\tENSG2.1\tTX2.1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    symbol_out = tmp_path / "symbol.bed"
    ensid_out = tmp_path / "ensid.bed"

    db.write_deduplong(gene_bed, symbol_out)
    db.write_deduplong(gene_bed, ensid_out, gene_key="ensid")

    assert len(symbol_out.read_text(encoding="utf-8").splitlines()) == 1
    assert len(ensid_out.read_text(encoding="utf-8").splitlines()) == 2


def test_write_tss_flank_region_unions(tmp_path):
    gene_bed = tmp_path / "genes.bed"
    gene_bed.write_text(
        "\n".join(
            [
                "chr1\t100\t200\tGeneA\t100\t+\tGENEA\tTX1",
                "chr1\t150\t260\tGeneA\t110\t+\tGENEA\tTX2",
                "chr1\t400\t500\tGeneB\t100\t-\tGENEB\tTX3",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    outputs = db.write_tss_flank_region_unions(
        gene_bed,
        tmp_path,
        promoter_bp=20,
        distal_bp=50,
        prefix="20bp",
    )

    assert outputs["promoter.up"] == tmp_path / "20bp.promoter.up.bed"
    assert outputs["promoter.up"].read_text(encoding="utf-8") == (
        "chr1\t80\t100\n"
        "chr1\t130\t150\n"
        "chr1\t500\t520\n"
    )
    assert outputs["promoter.down"].read_text(encoding="utf-8") == (
        "chr1\t101\t121\n"
        "chr1\t151\t171\n"
        "chr1\t479\t499\n"
    )
    assert outputs["dis5"].read_text(encoding="utf-8") == (
        "chr1\t50\t80\n"
        "chr1\t100\t130\n"
        "chr1\t520\t550\n"
    )
    assert outputs["dis3"].read_text(encoding="utf-8") == (
        "chr1\t220\t250\n"
        "chr1\t280\t310\n"
        "chr1\t350\t380\n"
    )


def _write_isoform_selection_gene_bed(tmp_path):
    gene_bed = tmp_path / "all.gene.bed"
    gene_bed.write_text(
        "\n".join(
            [
                "chr1\t100\t200\tGeneA\t100\t+\tGENEA.1\tTXA1.1\tprotein_coding\tGENEA\tTXA1",
                "chr1\t140\t300\tGeneA\t160\t+\tGENEA.1\tTXA2.1\tprotein_coding\tGENEA\tTXA2",
                "chr1\t400\t500\tGeneB\t100\t-\tGENEB.1\tTXB1.1\tprotein_coding\tGENEB\tTXB1",
                "chr1\t420\t650\tGeneB\t230\t-\tGENEB.1\tTXB2.1\tprotein_coding\tGENEB\tTXB2",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return gene_bed


def test_dedup_gencode_bed_isoid_falls_back_to_longest(tmp_path):
    gene_bed = _write_isoform_selection_gene_bed(tmp_path)
    ids = tmp_path / "ids.txt"
    ids.write_text("TXA1\n", encoding="utf-8")

    outputs = db.dedup_gencode_bed(
        "isoID",
        ids,
        tmp_path,
        gene_bed=gene_bed,
    )

    assert outputs["gene"] == tmp_path / "all.gene.dedupisoID.gene.bed"
    rows = outputs["gene"].read_text(encoding="utf-8").splitlines()
    assert [row.split("\t")[7] for row in rows] == ["TXA1.1", "TXB2.1"]
    assert outputs["tss"].exists()
    assert outputs["tes"].exists()


def test_dedup_gencode_bed_long_selectors(tmp_path):
    gene_bed = _write_isoform_selection_gene_bed(tmp_path)

    col5_outputs = db.dedup_gencode_bed(
        "longcol5", output_dir=tmp_path / "col5", gene_bed=gene_bed
    )
    interval_outputs = db.dedup_gencode_bed(
        "long", output_dir=tmp_path / "interval", gene_bed=gene_bed
    )
    default_outputs = db.dedup_gencode_bed(
        output_dir=tmp_path / "default", gene_bed=gene_bed
    )

    assert [
        row.split("\t")[7]
        for row in col5_outputs["gene"].read_text(encoding="utf-8").splitlines()
    ] == ["TXA2.1", "TXB2.1"]
    assert [
        row.split("\t")[7]
        for row in interval_outputs["gene"].read_text(encoding="utf-8").splitlines()
    ] == ["TXA2.1", "TXB2.1"]
    assert default_outputs["gene"].name == "all.deduplongcol5.gene.bed"


def test_parse_bp_accepts_m_and_mb():
    assert dedup._parse_bp("2m") == 2000000
    assert dedup._parse_bp("1.5mb") == 1500000
    assert regions._parse_bp("3m") == 3000000
    assert regions._parse_bp("2mb") == 2000000


def test_ensembl_gtf_url_uses_release_and_reference():
    assert ensembl_gtf_url("human", "100") == (
        "https://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/"
        "Homo_sapiens.GRCh38.100.gtf.gz"
    )
    assert ensembl_gtf_url("arabidopsis", "50") == (
        "https://ftp.ensemblgenomes.ebi.ac.uk/pub/release-50/plants/gtf/"
        "arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.50.gtf.gz"
    )


def test_ensembl_gtf_url_uses_cached_genomes_species(monkeypatch, tmp_path):
    monkeypatch.setenv("SJCAB_PEAK2ANNO_DB_PATH", str(tmp_path))
    cache = tmp_path / "ensembl" / "63" / "species.txt"
    cache.parent.mkdir(parents=True)
    cache.write_text(
        "species\tdivision\tassembly\n"
        "bigelowiella_natans\tEnsemblProtists\tBigna1\n",
        encoding="utf-8",
    )
    assert ensembl_gtf_url("bigelowiella_natans", "63") == (
        "https://ftp.ensemblgenomes.ebi.ac.uk/pub/release-63/protists/gtf/"
        "bigelowiella_natans/Bigelowiella_natans.Bigna1.63.gtf.gz"
    )


def test_filter_gencode_bed_isoid_omits_unmatched_genes(tmp_path):
    gene_bed = _write_isoform_selection_gene_bed(tmp_path)
    ids = tmp_path / "ids.txt"
    ids.write_text("TXA1\n", encoding="utf-8")

    outputs = db.filter_gencode_bed(
        "isoID",
        ids,
        tmp_path,
        gene_bed=gene_bed,
    )

    rows = outputs["gene"].read_text(encoding="utf-8").splitlines()
    assert [row.split("\t")[7] for row in rows] == ["TXA1.1"]


def test_dedup_gencode_bed_symbol_default_can_select_ensid_grouping(tmp_path):
    gene_bed = tmp_path / "my.bed"
    gene_bed.write_text(
        "\n".join(
            [
                "chr1\t100\t200\tSame\t100\t+\tENSG1.1\tTX1.1\tprotein_coding\tENSG1\tTX1",
                "chr1\t120\t240\tSame\t120\t+\tENSG2.1\tTX2.1\tprotein_coding\tENSG2\tTX2",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    ids = tmp_path / "ids.txt"
    ids.write_text("TX1\n", encoding="utf-8")

    symbol_outputs = db.dedup_gencode_bed(
        "isoID",
        ids,
        tmp_path / "symbol",
        gene_bed=gene_bed,
    )
    ensid_outputs = db.dedup_gencode_bed(
        "isoID",
        ids,
        tmp_path / "ensid",
        gene_bed=gene_bed,
        gene_key="ensid",
    )

    assert symbol_outputs["gene"].name == "my.dedupisoID.gene.bed"
    assert ensid_outputs["gene"].name == "my.dedupisoID.gene.bed"
    assert len(symbol_outputs["gene"].read_text(encoding="utf-8").splitlines()) == 1
    assert len(ensid_outputs["gene"].read_text(encoding="utf-8").splitlines()) == 2


def test_dedup_gencode_bed_species_version_default_prefix(monkeypatch, tmp_path):
    gene_bed = _write_isoform_selection_gene_bed(tmp_path)
    ids = tmp_path / "ids.txt"
    ids.write_text("TXA1\n", encoding="utf-8")

    def fake_registry_path(*args, **kwargs):
        return gene_bed

    monkeypatch.setattr(dedup, "registry_path", fake_registry_path)

    outputs = db.dedup_gencode_bed(
        "isoID",
        ids,
        tmp_path,
        species="mm10",
        version="vM22",
    )

    assert outputs["gene"] == tmp_path / "mm10.vM22.dedupisoID.gene.bed"


def test_dedup_gencode_bed_isoexp_uses_highest_expression(tmp_path):
    gene_bed = _write_isoform_selection_gene_bed(tmp_path)
    expression = tmp_path / "expr.tsv"
    expression.write_text(
        "isoform\texpression\n"
        "TXA1\t1\n"
        "TXA2\t9\n"
        "TXB1\t3\n",
        encoding="utf-8",
    )

    outputs = db.dedup_gencode_bed(
        "isoexp",
        expression,
        tmp_path,
        gene_bed=gene_bed,
    )

    rows = outputs["gene"].read_text(encoding="utf-8").splitlines()
    assert [row.split("\t")[7] for row in rows] == ["TXA2.1", "TXB1.1"]


def test_dedup_gencode_bed_peak_uses_max_promoter_peak_score(tmp_path):
    gene_bed = _write_isoform_selection_gene_bed(tmp_path)
    peaks = tmp_path / "peaks.bed"
    peaks.write_text(
        "chr1\t90\t95\tpeakA\t10\n"
        "chr1\t145\t150\tpeakB\t5\n",
        encoding="utf-8",
    )

    outputs = db.dedup_gencode_bed(
        "peak",
        peaks,
        tmp_path,
        gene_bed=gene_bed,
        promoter_bp=20,
    )

    rows = outputs["gene"].read_text(encoding="utf-8").splitlines()
    assert [row.split("\t")[7] for row in rows] == ["TXA1.1", "TXB2.1"]


def test_dedup_gencode_bed_peak_accepts_region_text(tmp_path):
    gene_bed = _write_isoform_selection_gene_bed(tmp_path)
    peaks = tmp_path / "peaks.txt"
    peaks.write_text(
        "region score\n"
        "chr1:90-95 10\n"
        "chr1_145_150 5\n",
        encoding="utf-8",
    )

    outputs = db.dedup_gencode_bed(
        "peak", peaks, tmp_path / "peak-text", gene_bed=gene_bed, promoter_bp=20
    )

    rows = outputs["gene"].read_text(encoding="utf-8").splitlines()
    assert [row.split("\t")[7] for row in rows] == ["TXA1.1", "TXB2.1"]


def test_filter_gencode_bed_perover_uses_promoter_overlap_percentage(tmp_path):
    gene_bed = _write_isoform_selection_gene_bed(tmp_path)
    active = tmp_path / "active.bed"
    active.write_text(
        "chr1\t90\t100\n"
        "chr1\t145\t146\n",
        encoding="utf-8",
    )

    outputs = db.filter_gencode_bed(
        "perover",
        active,
        tmp_path,
        gene_bed=gene_bed,
        promoter_bp=20,
    )

    rows = outputs["gene"].read_text(encoding="utf-8").splitlines()
    assert [row.split("\t")[7] for row in rows] == ["TXA1.1"]


def test_filter_gencode_bed_perover_accepts_region_text_without_header(tmp_path):
    gene_bed = _write_isoform_selection_gene_bed(tmp_path)
    active = tmp_path / "active.txt"
    active.write_text("chr1:90-100\nchr1/145/146\n", encoding="utf-8")

    outputs = db.filter_gencode_bed(
        "perover", active, tmp_path / "perover-text", gene_bed=gene_bed, promoter_bp=20
    )

    rows = outputs["gene"].read_text(encoding="utf-8").splitlines()
    assert [row.split("\t")[7] for row in rows] == ["TXA1.1"]


def test_dedup_gencode_bed_cli_accepts_gene_bed(monkeypatch, tmp_path):
    called = {}
    selector = tmp_path / "ids.txt"
    selector.write_text("TXA1\n", encoding="utf-8")
    gene_bed = tmp_path / "all.gene.bed"
    gene_bed.write_text("", encoding="utf-8")

    def fake_dedup(method, selector_path, **kwargs):
        called["method"] = method
        called["selector"] = selector_path
        called["kwargs"] = kwargs
        return {"gene": tmp_path / "dedupisoID.gene.bed"}

    monkeypatch.setattr(cli, "dedup_gencode_bed", fake_dedup)

    assert (
        cli.main(
            [
                "dedup-bed",
                "-b",
                str(gene_bed),
                "-m",
                "isoID",
                "-i",
                str(selector),
                "-o",
                str(tmp_path),
                "--exclusive",
            ]
        )
        == 0
    )
    assert called == {
        "method": "isoID",
        "selector": str(selector),
        "kwargs": {
            "output_dir": str(tmp_path),
            "gene_bed": str(gene_bed),
            "species": None,
            "version": "default",
            "data_dir": None,
            "promoter_bp": "2kb",
            "inclusive": False,
            "output_prefix": None,
            "gene_key": "symbol",
        },
    }
