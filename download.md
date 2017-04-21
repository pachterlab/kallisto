---
layout: page
title: "Download"
group: navigation
---

{% include JB/setup %}

If using a Mac the easiest way to download and install __kallisto__ is via brew with the commands:

```
ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
brew tap homebrew/science
brew install kallisto
```

For other operating systems download executables or source from the links below and follow the provided instructions.

#### Releases

The __kallisto__ GitHub repository is [here](http://github.com/pachterlab/kallisto).

<table class="table">
  <thead>
    <tr>
      <th style="text-align: left">Version</th>
      <th>Date</th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>

{% for post in site.categories.releases %}
    <tr>
    	<td>Release notes: <a href="{{ site.url }}/kallisto/{{ post.url }}">{{ post.version }}</a></td>
    	<td><span class="entry-date"><time datetime="{{ post.date | date_to_xmlschema }}">{{ post.date | date: "%B %d, %Y" }}</time></span></td>

        <td><a href="https://github.com/pachterlab/kallisto/releases/download/{{ post.version }}/kallisto_mac-{{ post.version }}.tar.gz">Mac</a></td>
        <td>{% if post.version != 'v0.42.2' %} <a href="https://github.com/pachterlab/kallisto/releases/download/{{ post.version }}/kallisto_linux-{{ post.version }}.tar.gz">Linux</a> {% endif %} </td>
        <td>{% if post.version >= 'v0.43.0' %} <a href="https://github.com/pachterlab/kallisto/releases/download/{{ post.version }}/kallisto_windows-{{ post.version }}.zip">Windows</a> {% endif %} </td>
        <td><a href="https://github.com/pachterlab/kallisto/archive/{{ post.version }}.tar.gz">Source</a></td>
    </tr>
{% endfor %}
</table>

#### Transcriptomes

Commonly used transcriptome FASTA files are available [here](http://bio.math.berkeley.edu/kallisto/transcriptomes/).

#### Licence

Kallisto is distributed under a non-commercial license reproduced below. For commercial use, please contact Terri Sale at the Office of Technology Licensing, UC Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley, CA 94720-1620, (510) 643-4219, terri.sale@berkeley.edu.


~~~
Copyright ©2015. The Regents of the University of California (Regents). All Rights Reserved. Permission to use, copy, modify, and distribute this software and its documentation for educational and research not-for-profit purposes, without fee and without a signed licensing agreement, is hereby granted, provided that the above copyright notice, this paragraph and the following two paragraphs appear in all copies, modifications, and distributions. Contact The Office of Technology Licensing, UC Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley, CA 94720-1620, (510) 643-7201, for commercial licensing opportunities.

Created by Nicolas Bray, Harold Pimentel, Pall Melsted and Lior Pachter, University of California, Berkeley

IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
~~~
