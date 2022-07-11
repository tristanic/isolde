$(document).ready(function () {
  $('a[href^="http://"], a[href^="https://"], a[href^="help:"]').not('a[class*=internal]').attr('target', '_blank');
});
