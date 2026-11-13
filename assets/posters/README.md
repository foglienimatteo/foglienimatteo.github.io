Posters and other PDFs shown as a preview card on /Publications/.

To add one:

1. drop the PDF here, e.g. `ISC26_Poster.pdf`
2. `bin/pdf-preview assets/posters/ISC26_Poster.pdf`
   (writes `ISC26_Poster.jpg` next to it)
3. reference both from `Publications/index.qmd`:

       {{< pdf-embed pdf="/assets/posters/ISC26_Poster.pdf"
                     image="/assets/posters/ISC26_Poster.jpg"
                     alt="..."
                     caption="..." >}}

Size and sharpness are two separate knobs, in two different places:

  - how big the card is drawn: `width` / `height` on the shortcode (CSS
    lengths). Without either, the card is capped at 34rem high.
  - how sharp the image is: `-w` (raster pixels) and `-q` (JPEG quality) on
    bin/pdf-preview, because those are properties of the file on disk.
    Render at roughly twice the `width` you draw at.

/CV/ uses the same shortcode at `width="36rem"`, from a raster made with
`bin/pdf-preview -w 1400 -q 85`.

Poster PDFs exported from Illustrator or Inkscape are often 50-200MB. Git
keeps every binary forever, so shrink it first, e.g.

    gs -sDEVICE=pdfwrite -dPDFSETTINGS=/ebook -dNOPAUSE -dQUIET -dBATCH \
       -sOutputFile=small.pdf big.pdf
