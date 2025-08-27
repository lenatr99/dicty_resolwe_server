/**
 * Return deserialized response (type T) or null if response is not ok.
 * @param response Deserialized response json (as T).
 */
export const deserializeResponse = async <T>(response: Response): Promise<T> => {
    const responseJson = await response.json();
    return (responseJson.results != null ? responseJson.results : responseJson) as T;
};
