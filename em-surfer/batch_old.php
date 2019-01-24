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
								<td class="text6"><p>&nbsp;</p>Submit your EMDB id list in a batch mode<p>&nbsp;</p></td>
							</tr>
							
							<tr>
								<td>
									<form action="cgi-bin/searchbatch.php" method="post"
										enctype="multipart/form-data">

										<script><!--
		function qs(el) {if (window.RegExp && window.encodeURIComponent) {var qe=encodeURIComponent(document.f.q.value);if (el.href.indexOf("q=")!=-1) {el.href=el.href.replace(new RegExp("q=[^&$]*"),"q="+qe);} else {el.href+="&q="+qe;}}return 1;}
		// -->
		  </script>
										<table width="748">
											<td width="461"></c:if>

												<table width="220">
													<tr>
														<td width="220" class="text4"><p>&nbsp;</p>
															<p>Step 1 (Representation)</p></td>
													</tr>
												</table>
												<table width="650" border="1">
													<tr>
														<td width="215">&nbsp;&nbsp;&nbsp;&nbsp;Contour shape representation:</td>
														<td width="419"><select name="rep" id="rep">
														<option value="recommend">EMDB contour</option>
														<option value="recommendonethird">EMDB contour + 1/3 core</option>
														<option value="recommendtwothird" selected="selected">EMDB contour + 2/3 core</option>
														<option value="recommendonethirdtwothird">EMDB contour + 1/3 + 2/3 core</option>
												</select></td>
													</tr>
												</table>
												
												<table width="220">
													<tr>
														<td width="196" class="text4"><p>&nbsp;</p>
															<p>Step 2 (Query EMDB ID list)</p></td>
													</tr>
												</table>


												<table width="650" border="1">
													<tr>
														<td>&nbsp;&nbsp;&nbsp;&nbsp;ID list:</td>
														<td><textarea name="listarea" id="listarea" cols="40"
																rows="7">1010&#13;1884&#13;5502</textarea></td>
														<td><strong>e.g.</strong><br> <br>1010<br>1884<br>5502<br></td>
													</tr>

													<tr>
														<td colspan="3" align="center"><strong>Or</strong>
														</td>
													</tr>

													<tr>
														<td>&nbsp;&nbsp;&nbsp;&nbsp;Upload a file of the EMDB id list:</td>
														<td colspan="2"><input type="file" name="listfile"
															id="listfile" />
															
															<a href="idlist.dat">An example ID list file.</a>
															</td>
													</tr>
													
													
												</table>
												

												<table>
													<tr>
														<td class="text4"><p>&nbsp;</p>
															<p>Step 3 (Filter)</p></td>
													</tr>
												</table>
												<table width="650" border="1">
													<tr>
													
														<td>&nbsp;<b>Volume filter:</b>
												<br>
												(The volume of the EM entry in the database
												<br> 
												is between 0.8 and 1.2 times the volume of
												<br>
												the query entry if ON, or else if OFF)
												<br>
												</td>
														<td><input type="radio" name="volumeFilter"
															id="volumeFilter" value="on" checked="checked"/>
															ON&nbsp;&nbsp;&nbsp;&nbsp; <input type="radio"
															name="volumeFilter" id="volumeFilter" value="off"/> OFF
														</td>
													</tr>
												</table>

												<table width="135">
													<tr>
														<td width="127" class="text4"><p>&nbsp;</p>
															<p>Step 4 (Result)</p></td>
													</tr>
												</table>
												<table width="650" border="1">
													<tr>
														<td width="215">&nbsp;&nbsp;&nbsp;&nbsp;Top results shown:</td>
														<td width="419"><select name="num" id="num">
																<option value="10">10</option>
																<option value="20" selected="selected">20</option>
																<option value="30">30</option>
														</select></td>
													</tr>
												</table>

												<p>&nbsp;</p>
												<table width="650">

													<tr>
														<td align="center"><input type="submit" value="Submit" />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<input
															type="reset" value="Reset" /> <p>&nbsp;</p></td>
													</tr>
												</table>
										
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
