class Nesoni < Formula
  homepage "http://www.bioinformatics.net.au/software.nesoni.shtml"
  url "https://github.com/Victorian-Bioinformatics-Consortium/nesoni/archive/v0.130.tar.gz"
  version "0.130"
  sha1 "4077a992dfc519f3bd6a73fe76e4b90eb5ccc705"

  # For matplotlib
  depends_on "freetype"
  depends_on "samtools"
  depends_on "shrimp"

  resource "numpy" do
    url "https://pypi.python.org/packages/source/n/numpy/numpy-1.9.2.tar.gz"
    sha1 "86b4414cd01c4244141c59ea476ca8fdad8e9be2"
  end

  resource "matplotlib" do
    url "https://pypi.python.org/packages/source/m/matplotlib/matplotlib-1.4.3.tar.gz"
    sha1 "8b24472780a23e686135dd08ec1bc5ef88db1979"
  end

  resource "biopython" do
    url "https://pypi.python.org/packages/source/b/biopython/biopython-1.65.tar.gz"
    sha1 "16dc84622a707118b3cd48c09807810de61d2321"
  end


  def install
    ENV.prepend_create_path "PYTHONPATH", libexec/"vendor/lib/python2.7/site-packages"
    %w[numpy matplotlib biopython].each do |r|
      resource(r).stage do
        system "python", *Language::Python.setup_install_args(libexec/"vendor")
      end
    end

    ENV.prepend_create_path "PYTHONPATH", libexec/"lib/python2.7/site-packages"
    system "python", *Language::Python.setup_install_args(libexec)

    bin.install Dir["#{libexec}/bin/*"]
    bin.env_script_all_files(libexec/"bin", :PYTHONPATH => ENV["PYTHONPATH"])
  end


  def caveats
    <<-EOS.undent
      
      1) Nesoni has been installed without the dependencies required for VCF calling:
        * Picard
        * Freebayes
        * SplitsTree4
      2) We also ignore R libraries (we don't do RNA-Seq)
    EOS
  end


  test do
    system "#{bin}/nesoni", "clip"
  end
end
