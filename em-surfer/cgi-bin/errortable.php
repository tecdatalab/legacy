<?php

// Defines a function that can be called to generate a table-like formatted
// error, that complies with roughly the same format as other pages
function generate_error_table($error_title, $error_message) {
echo <<<END
<table width="100%" border="0" cellspacing="0" cellpadding="0">
	<tr>
		<td>
			<table width="700" border="0" align="center" cellpadding="0"
				cellspacing="0">
				<tr valign="top">
					<td width="700">
						<table width="100%" border="0" cellspacing="0" cellpadding="0">
							<tr>
								<td><img src="images/spacer.gif" alt="" width="100" height="35" />
								</td>
							</tr>
							<tr>
								<td class="text4">$error_title</td>
							</tr>
							<tr>
								<td>&nbsp;</td>
							</tr>
							<tr>
								<td>

									<table width="610" border="0" align="center" cellpadding="0"
										cellspacing="0">
										<tbody>
											<tr>
												<td colspan="4"><table cellspacing="0" cellpadding="0">
														<tr>
															<td width="610">
															$error_message
															<p>&nbsp;</p>
															<p>&nbsp;</p></td>
														</tr>
														<br />
													</table></td>
											</tr>
										</tbody>
									</table>
								</td>
							</tr>
							<tr>
								<td>&nbsp;</td>
							</tr>
						</table>
					</td>
				</tr>
			</table>
		</td>
	</tr>
</table>
END;
}
?>

