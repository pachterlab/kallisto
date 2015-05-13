---
layout: page
title: "Download"
group: navigation
---

{% include JB/setup %}

#### Source

The __kallisto__ source code is available
[here](http://github.com/pachterlab/kallisto).

#### Releases

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
    	<td><a href="{{ site.url }}{{ post.url }}">{{ post.version }}</a></td>
    	<td><span class="entry-date"><time datetime="{{ post.date | date_to_xmlschema }}">{{ post.date | date: "%B %d, %Y" }}</time></span></td>

        <td><a href="https://github.com/pachterlab/kallisto/releases/download/{{ post.version }}/kallisto_mac-{{ post.version }}.tar.gz">Mac</a></td>
        <td><a href="https://github.com/pachterlab/kallisto/releases/download/{{ post.version }}/kallisto_linux-{{ post.version }}.tar.gz">Linux</a></td>
        <td><a href="https://github.com/pachterlab/kallisto/archive/{{ post.version }}.tar.gz">Source</a></td>
    </tr>
{% endfor %} 
</table>


#### Transcriptomes

Commonly used transcriptome FASTA files are available [here](http://bio.math.berkeley.edu/kallisto/transcriptomes/).

### Licence

Kallisto is distributed under a non-commercial [license](license.html). For commercial use, contact The Office of Technology Licensing, UC Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley, CA 94720-1620, (510) 643-7201
