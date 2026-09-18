Posters and other PDFs shown as a preview card on /Publications/.

To add one:

1. drop the PDF here, e.g. `ISC26_Poster.pdf`
2. `bin/pdf-preview assets/posters/ISC26_Poster.pdf`
   (writes `ISC26_Poster.jpg` next to it)
3. reference both from `_pages/Work.md`:

       {% include pdf-embed.html
            pdf="/assets/posters/ISC26_Poster.pdf"
            image="/assets/posters/ISC26_Poster.jpg"
            alt="..."
            caption="..." %}

The include renders nothing until the PDF actually exists, so an unfinished
entry never shows a broken image on the live site.

Poster PDFs exported from Illustrator or Inkscape are often 50-200MB. Git
keeps every binary forever, so shrink it first, e.g.

    gs -sDEVICE=pdfwrite -dPDFSETTINGS=/ebook -dNOPAUSE -dQUIET -dBATCH \
       -sOutputFile=small.pdf big.pdf
