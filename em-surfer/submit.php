<?php
include ('header.php');
?>

<script
	type='text/javascript' src='jquery.all.js'></script>
<link
	rel="stylesheet" type="text/css" href="jquery.autocomplete.css" />
<script type="text/javascript">
$().ready(function() {
        $("#q").autocomplete("searchPDB.php", {
                width: 200,
                selectFirst: false
        });
});

var validateForm;
validateForm = function() {
	 if ( form.q.value == "" && form.pdbfile.value == "") {
	  	alert( "\n Please input a valid protein structure ID or upload your own structure file.");
	  	form.q.focus();
	 	return false; // Stop the form processing.
	 }

	 //if ( form.q.value != "" && form.pdbfile.value != "") {
	 //	  	alert( "\n Please input an EMDB ID or upload your own map/mrc file, not both.");
	 //  	form.q.focus();
	//	 	return false; // Stop the form processing.
	// }

	if ( form.q.value != "" && form.pdbfile.value == "") {
		if ( !/^[0-9]{4,4}$/.test(form.q.value )) {
			  alert( "\n Invalid ID format or length. Input protein id's must be numbers(0-9) and 4 numbers long" );
			  form.q.focus();
			 return false; // Stop the form processing.
		}
	}
	 
	return true; // Submit the form if your codintion is achieved.
}
</script>

<table width="100%" align="center" border="0" cellspacing="0" cellpadding="0">
	<tr>
		<td>
			<table width="760" border="0" align="center" cellpadding="0"
				cellspacing="0">
				<tr valign="top">
					<td width="760" align="center">
						<table width="100%" border="0">
							<tr>
								<td><p>&nbsp;</p>
									<p class="text6">Submit an EM map</p> Please refer to the tips
									below when uploading a file:
									<ul>
										<li>It is possible that the EMDB ID already exists in the
											database. Try to use the search box first. </li>
										<li>If you benchmark our program, please
										<u><a href="batch.php">access the page
										</a> 
										</u>
										.
										</li>
									</ul>
									<p class="text4"> (For the purpose of review, please directly click the <em>Submit</em> button to get the search result, 
									<br>
									with the default options provided.) </p>
								</td>
							</tr>

							<tr>
								<td>
									<form method="post" action="cgi-bin/search.php"
										enctype="multipart/form-data"
										onsubmit="return validateForm()">

										<script><!--
		function qs(el) {if (window.RegExp && window.encodeURIComponent) {var qe=encodeURIComponent(document.f.q.value);if (el.href.indexOf("q=")!=-1) {el.href=el.href.replace(new RegExp("q=[^&$]*"),"q="+qe);} else {el.href+="&q="+qe;}}return 1;}
		// -->
		  </script>
										<p class="text4">Step 1 (Representation)</p>
										<table width="680" border="1" align="center">
											<tr>
												<td width="240">&nbsp;&nbsp;&nbsp;&nbsp;Contour shape representation:</td>
												<td width="440"><select name="rep" id="rep">
														<option value="recommend" selected="selected">EMDB contour</option>
														<option value="recommendonethird">EMDB contour + 1/3 core</option>
														<option value="recommendtwothird">EMDB contour + 2/3 core</option>
														<option value="recommendonethirdtwothird">EMDB contour + 1/3 + 2/3 core</option>
														<option value="recommendonestd">EMDB contour + 1 std dev</option>
												</select></td>
											</tr>
										</table>
										
										<p>&nbsp;</p><p class="text4">Step 2 (Query entry)</p>
										<table width="680" border="1" align="center">
											<tr>
												<td width="220">
													<p><strong> Enter 4-digit EMDB entry ID: <br> e.g.,</strong> ID: 1884<br></p>
													<p>
														<input name="q" type="text" value="1884" id="q" maxlength="4" />
													</p> 
												</td>
												<td width="40" align="center"><strong>Or</strong>
												</td>
												<td width="420">
					`								<p>Upload an EM map (.map or .mrc) file
													(<a href="http://kiharalab.org/em-surfer/tutorial.php#upload_troubleshooting">Upload troubleshooting</a>):</p>
													<p>
														<input type="file" name="pdbfile" id="pdbfile" />
												        <a href="emd_1884.map">An example file.</a>
													</p>
													<p class="text4">and</p>
													<strong>Recommended contour level:</strong>
													<input name="contour" type="text" value="3.16" id="contour" />
									                <strong>e.g.,</strong> 3.16
												</td>
											</tr>
										</table>

										<p>&nbsp;</p><p class="text4">Step 3 (Filter)</p>
										<table width="680" border="1" align="center">
											<tr>
												<td width="420">&nbsp;<b>Volume filter:</b>
												<br>
												(The volume of the EM entry in the database is between 0.8 and 1.2
												<br>
												times the volume of the query entry if ON, or else if OFF)
												<br>
												</td>
												<td width="260"><input type="radio" name="volumeFilter"
													id="volumeFilter" value="on" checked="checked"/>
													ON&nbsp;&nbsp;&nbsp;&nbsp; <input type="radio"
													name="volumeFilter" id="volumeFilter" value="off"/> OFF</td>
											</tr>
											<tr>
												<td width="420">&nbsp;<b>Resolution filter:</b>
												<br>
												(The query is only compared against maps in this resolution range.
												<br> 
												If both are left blank, no filtering is applied. If only one is
												<br>
												provided, it will be the only restriction imposed)
												<br>
												</td>
												<td width="260">
												Min: 
												<input type="number" name="minResolutionFilter" id="minResolutionFilter"
													min="0" class="numericinput" step="0.0001"/>
												Max:
												<input type="number" name="maxResolutionFilter" id="maxResolutionFilter"
													min="0" class="numericinput" step="0.0001"/>
												</td>
											</tr>
										</table>
										<p>&nbsp;</p>
										<table width="680" align="center">
											<tr>
												<td align="center"><input type="submit" value="Submit" />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<input
													type="reset" value="Reset" />
													<p>&nbsp;</p></td>
											</tr>
										</table>
										
									</form>
								</td>
							</tr>
						</table>
					</td>
				</tr>
			</table>
		</td>
	</tr>
</table>
<?php
include ('tailer.php');
?>
